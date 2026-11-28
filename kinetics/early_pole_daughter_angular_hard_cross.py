"""Four-corner contrast in hard depth and scalar angular sampling."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from early_pole_daughter_angular_held_chart_compare import compare


def cross(a,b):
    if a['angular_mesh_sha256']!=b['angular_mesh_sha256'] or a['exposure']!=b['exposure']:
        raise ValueError('Hard-depth runs do not share scalar sampling or reporting times.')
    if any(set(v['paired_count_histories'])!={'angular','ordinary'} for v in (a,b)):
        raise ValueError('Unexpected corner receipt.')
    t=np.array(a['exposure'])
    x,y,z,u=(np.asarray(d['paired_count_histories'][k]) for d in (a,b) for k in ('ordinary','angular'))
    H=y-x
    R=u-z
    groups=np.array(b['angular_grouped_signed_gains'])-np.array(a['angular_grouped_signed_gains'])
    if groups.shape!=(len(t),8) or np.max(abs(groups.sum(axis=1)-(u-y)))>4e-17:
        raise ValueError('Regrouped hard contrast failed scalar chart closure.')
    names=('hard1_cold32','hard1_angular78','hard4_cold32','hard4_angular78')
    sources=dict(zip(names,(x,y,z,u),strict=True))
    shifts={k:float(v[-1]/t[-1]) for k,v in sources.items()}
    if any(not np.isfinite(v) for v in shifts.values()):
        raise ValueError('Invalid corner histories.')
    crossings={k:[t[i:i+2].tolist() for i in np.flatnonzero((v[:-1]<0)&(v[1:]>=0))]
               for k,v in sources.items()}
    return dict(scope="Actual held full local intervention at two hard depths and two scalar samplings, no continuum bound.",
                exposure=t.tolist(),angular_mesh_sha256=a['angular_mesh_sha256'],
                four_paired_means=shifts,negative_to_nonnegative_brackets=crossings,
                hard1_angular_shift=float(H[-1]/t[-1]),hard4_angular_shift=float(R[-1]/t[-1]),
                cold32_hard_shift=float((z[-1]-x[-1])/t[-1]),
                angular78_hard_shift=float((u[-1]-y[-1])/t[-1]),
                mixed_four_corner_correction=float((R[-1]-H[-1])/t[-1]),
                mixed_four_corner_history=(R-H).tolist(),
                hard_shift_in_angular_grouped_counts=groups.tolist(),
                hard1_receipt=a,hard4_receipt=b)


def read(p):
    if p.endswith('.gz'):
        with gzip.open(p,'rt') as f:
            return json.load(f)
    return json.loads(Path(p).read_text())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('hard1-angular','hard1-ordinary','hard4-angular','hard4-ordinary'):
        p.add_argument('--'+k,nargs=2,required=True)
    p.add_argument('--initializer',required=True)
    a=p.parse_args()
    fixtures={d['depth']:d for d in read(a.initializer)}
    args=(a.hard1_angular,a.hard1_ordinary,a.hard4_angular,a.hard4_ordinary)
    values=[[read(v) for v in q] for q in args]
    for d,i in ((1,0),(4,2)):
        if any(v.get('aperture_subdivision',1)!=d for q in values[i:i+2] for v in q):
            raise ValueError('Unexpected hard benchmark depth.')
    receipt=cross(compare(values[0],values[1],fixtures[1]),
                  compare(values[2],values[3],fixtures[4]))
    early,late=values[0],values[2]
    np.testing.assert_array_equal(early[0]['cold_momentum_nodes_over_T'],late[0]['cold_momentum_nodes_over_T'])
    np.testing.assert_array_equal(early[0]['cold_mode_count_history'][0],late[0]['cold_mode_count_history'][0])
    a,b=(np.asarray(d[1]['cold_mode_count_history'][-1])-np.asarray(d[0]['cold_mode_count_history'][-1])
         for d in (early,late))
    local=late[1]['local_archive_support']
    if len(local['indices'])!=8 or len(set(local['indices']))!=8:
        raise ValueError('Hard-4 local intervention did not use all eight nested daughter cells.')
    receipt['angular_terminal_modes_smaller_at_hard4']=int(np.count_nonzero(b<a))
    receipt['angular_terminal_mode_gains']={'hard1':a.tolist(),'hard4':b.tolist()}
    receipt['hard4_local_support']=local
    receipt['energy_defect_max']=max(abs(z['energy_defect']) for q in values for z in q)
    receipt['entropy_defect_max']=max(abs(z['entropy_defect']) for q in values for z in q)
    print(json.dumps(receipt,indent=2,allow_nan=False))
