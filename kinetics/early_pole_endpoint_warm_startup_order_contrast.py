"""Compare actual warm quadrature startup endpoints on the same archived clock."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_chain_receipt import run as gaussian_audit
from early_pole_endpoint_gauss_warm_comparison import gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load
from early_pole_endpoint_targeted_warm_history import audit as targeted_audit


def run(native,native_chain,reports,states):
    x=.0075
    ref,sha=read(native)
    card=native_chain['validated_chain'][0]
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    values=np.asarray(ref['endpoint_occupation'])
    if (card['uncompressed_state_sha256']!=sha or (ref['grid'],ref['charge_div'],ref['scalar_div'],ref['cut'],
        ref['x_end'],ref['steps'],ref['counts'])!=('fine',2,4,1e-5,x,4096,[d.ns,d.nf,d.np])
        or values.shape!=(d.n,) or np.min(values)<0. or not np.all(np.isfinite(values)) or np.max(values[d.ns:])>1.):
        raise ArithmeticError('Cannot identify the full native S4 startup clock and state')
    mass=d.state(x,values)[0]
    labels=('G2','G3','G4','G5','G6','targeted_2xG4')
    rows=[]
    for label,order,r,s in zip(labels,(None,3,4,5,6,None),reports,states,strict=True):
        a=targeted_audit([r],[s],prefix=True) if label.startswith('target') else gaussian_audit([r],[s],prefix=True,order=order)
        q,_=read(s)
        z,h,_=gauss(q,d,mass,x,expected_steps=4096)
        rows.append(dict(chart=label,source=a['checked_epoch_pairs'][0],fixed_mass_count_and_pressures=z.tolist(),
                         conditional_pressure_envelope=h))
    pairs=[]
    for i,j in ((0,1),(1,2),(2,3),(3,4),(4,5),(3,5),(2,4)):
        a,b=[np.asarray(rows[k]['fixed_mass_count_and_pressures']) for k in (i,j)]
        l,r=[rows[k]['conditional_pressure_envelope']['pressure_total_min_max'] for k in (i,j)]
        bounds=[r[0]/l[1]-1.,r[1]/l[0]-1.]
        if not bounds[0]-2e-14<=b[1]/a[1]-1.<=bounds[1]+2e-14:
            raise ArithmeticError('Reconstructed startup pressure escaped its joint conditional bracket')
        pairs.append(dict(left=labels[i],right=labels[j],right_over_left_minus_one=(b/a-1.).tolist(),
                          pressure_relative_conditional_bounds=bounds))
    return dict(x=x,steps=[2048,4096],fixed_mass=mass,native_state_sha256_uncompressed=sha,
                fields=['count','pressure','warm_pressure'],charts=rows,contrasts=pairs,
                scope='Startup only at a fixed mass; this verifies a prefix, not subsequent warm evolution or a physical pressure prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',type=Path,required=True)
    p.add_argument('--native-chain',type=Path,required=True)
    for k in ('reports','states'):
        p.add_argument('--'+k,type=Path,nargs=6,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.native,json.loads(a.native_chain.read_text()),a.reports,a.states),indent=2,allow_nan=False))
