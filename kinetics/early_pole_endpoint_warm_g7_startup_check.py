"""Check the actual G6 and G7 warm startup pressure at a common native mass."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_chain_receipt import run as chain
from early_pole_endpoint_gauss_warm_comparison import gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load


def run(native,receipt,reports,states):
    x=.0075
    c,h=read(native)
    q=receipt['validated_chain'][0]
    b=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    o=np.asarray(c['endpoint_occupation'])
    if h!=q['uncompressed_state_sha256'] or (c['x_end'],c['steps'])!=(x,4096) or o.shape!=(b.n,):
        raise ArithmeticError('The common S4 startup state has changed')
    m=b.state(x,o)[0]
    rows=[]
    raw=[]
    for k,f,s in zip((6,7),reports,states,strict=True):
        audit=chain([f],[s],prefix=True,order=k)
        state,_=read(s)
        v,e,_=gauss(state,b,m,x,expected_steps=4096)
        rows.append(dict(order=k,source=audit['checked_epoch_pairs'][0],charted_moments=v.tolist(),
                         conditioned_pressure_interval=e['pressure_total_min_max']))
        raw.append(v)
    a,b=[v['conditioned_pressure_interval'] for v in rows]
    bound=[b[0]/a[1]-1.,b[1]/a[0]-1.]
    g=raw[1]/raw[0]-1.
    if not bound[0]-2e-14<=g[1]<=bound[1]+2e-14:
        raise ArithmeticError('Measured pressure escaped its independent conditioned interval')
    return dict(x=x,steps=[2048,4096],fixed_mass=m,fields=['count','pressure','warm_pressure'],
                charts=rows,G7_over_G6_minus_one=g.tolist(),conditioned_pressure_contrast_bound=bound,
                scope='Two completed startup epochs only, not full warm history, numerical convergence or physical prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',type=Path,required=True)
    p.add_argument('--native-chain',type=Path,required=True)
    for k in ('reports','states'):
        p.add_argument('--'+k,nargs=2,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.native,json.loads(a.native_chain.read_text()),a.reports,a.states),indent=2,allow_nan=False))
