"""Audit interaction between finite time and cold momentum quadrature errors."""
import argparse
import json

import numpy as np

from early_pole_joined_temporal_refinement_reference import select
from early_pole_first_window_cold_bin_refinement_reference import compare as check_mesh


def solve(a,b,c,d):
    ca=check_mesh([a,c])
    cb=check_mesh([b,d])
    method=('grid','projection','duration','pair','method','strength','compton_strength','same_soft_strength',
            'conversion_strength','screening_update','response_smear','aperture_grid','hard_edge_sets','elastic_energy_budget')
    if any(a[k]!=b[k] for k in method):
        raise ValueError('Time refinement must leave the physical controls unchanged.')
    qn=lambda q:len(q['snapshots'][0]['cold_bins'])
    if qn(a)!=qn(b) or qn(c)!=qn(d) or qn(c)<=qn(a):
        raise ValueError('Same cold meshes must be compared at both widths.')
    if b['steps']%a['steps'] or b['steps']<=a['steps']:
        raise ValueError('The coarse time grid must strictly nest in the fine grid.')
    for x,y in ((a,b),(c,d)):
        r,s=(a['snapshots'][0] for a in (x,y))
        for k in ('hard_moments','cold_bins','instantaneous_source_by_q'):
            np.testing.assert_array_equal(r[k],s[k])
    ratio=b['steps']//a['steps']
    v=[]
    step=a['steps']
    for k,(x,skip) in enumerate(((a,1),(b,ratio),(c,1),(d,ratio))):
        s=np.asarray(x['state_history'])[::skip]
        assert len(s)==step+1
        v.append(s-s[0])
    A,B,C,D=v
    mesh=np.array([C-A,D-B])
    time=np.array([B-A,D-C])
    mixed=mesh[1]-mesh[0]
    assert np.array_equal(mixed,time[1]-time[0])
    final=np.array([A[-1],B[-1],C[-1],D[-1]])
    r=a['duration']
    return dict(duration=r,time_grids=[a['steps'],b['steps']],cold_grids=[qn(a),qn(c)],
                projection=a['projection'],labels=['coarse_time_coarse_cold','fine_time_coarse_cold',
                'coarse_time_fine_cold','fine_time_fine_cold'],
                means=(final/r).tolist(),
                fractional_cold_mesh_effect_by_time_grid=[float(mesh[0,-1]/A[-1]),float(mesh[1,-1]/B[-1])],
                fractional_time_effect_by_cold_mesh=[float(time[0,-1]/A[-1]),float(time[1,-1]/C[-1])],
                signed_mixed_gain=float(mixed[-1]),max_abs_mixed_gain=float(np.max(abs(mixed))),
                signed_mixed_fraction_of_coarsest_gain=float(mixed[-1]/A[-1]),
                cold_mesh_audits=[ca,cb],scope='Matched conditional first-window calculation; finite cross difference, not a continuum estimate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('paths',nargs=4)
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    a=p.parse_args()
    print(json.dumps(solve(*(select(k,a.projection) for k in a.paths)),indent=2,allow_nan=False))
