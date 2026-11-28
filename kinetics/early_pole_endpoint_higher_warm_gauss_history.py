"""Advance independent higher-order warm Gauss nodes on the held cold and charge grid."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_stencil_remap_pilot import X,load
from early_pole_endpoint_warm_positive_gauss_ladder import WarmGauss
from early_pole_resolved_parent_reference import S


def run(order,end,steps,restart=None):
    d=WarmGauss(load('fine')[0],order)
    metadata=dict(grid='fine',charge_div=2,scalar_div=4,cut=1e-5,gauss_floor=.01,reference_x=X,
                  nodes=[d.ns,d.nf,d.np],scalar_cell_count=d.old_count,scalar_nodal_rule='gauss-volume-v2-general',
                  warm_gauss_order=order)
    return run_chart(d,metadata,end,steps,restart,scope='Independent positive warm Gauss nodal continuation, not a piecewise-constant archive.')


def run_chart(d,metadata,end,steps,restart,*,scope):
    z=None
    x0=X
    previous=None
    if restart is not None:
        raw=restart.read_bytes()
        p=json.loads(raw)
        if any(p.get(k)!=v for k,v in metadata.items()):
            raise ArithmeticError('Restart must preserve every chart element')
        z=np.asarray(p['nodal_occupation'])
        if z.shape!=(d.n,) or not np.all(np.isfinite(z)) or min(z)<0 or max(z[d.ns:])>1:
            raise ArithmeticError('Invalid restart occupations')
        x0=p['x_end']
        previous=hashlib.sha256(raw).hexdigest()
    if not X<=x0<end or steps<8 or steps%8:
        raise ValueError('Invalid clock segments')
    if (x0,end,steps) not in ((X,.0075,2048),(.0075,.01,512),(.01,.02,1024),(.02,.04,1024),(.04,.05,1024)):
        raise ValueError('Use the matched warm clock with replayable stage starts')
    coarse,_=drive(d,end,steps,x_start=x0,initial=z,stepper=step)
    fine,y=drive(d,end,2*steps,x_start=x0,initial=z,stepper=step)
    r=abs(coarse[-1]['scalar_yield']/fine[-1]['scalar_yield']-1)
    count=np.add.reduceat(d.ds.w*y[:d.ns]/S,np.r_[0,np.cumsum(d.repeats)[:-1]])
    if (r>1e-3 or max(abs(a['continuity_residual']) for a in coarse+fine)>3e-10
            or abs(fsum(count)/fine[-1]['scalar_yield']-1)>4e-14):
        raise ArithmeticError('Warm Gauss history violates its continuity, time or number gate')
    output=dict(metadata,x_end=end,steps=2*steps,actual_previous_raw_sha256=previous,
                nodal_occupation=y.tolist())
    return dict(metadata,starting_x=x0,x_end=end,coarse_steps=steps,checkpoint_trajectories=[coarse,fine],
                relative_time_drift=r,scalar_cell_yields=count.tolist(),scalar_pole_mass=d.state(end,y)[0],
                actual_previous_raw_sha256=previous,
                scope=scope),output


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--warm-order',type=int,choices=(3,4,5,6,7,8),default=3)
    p.add_argument('--x-end',type=float,default=.0075)
    p.add_argument('--steps',type=int,default=2048)
    p.add_argument('--restart',type=Path)
    p.add_argument('--state',type=Path,required=True)
    a=p.parse_args()
    r,z=run(a.warm_order,a.x_end,a.steps,a.restart)
    temp=a.state.with_suffix('.tmp')
    temp.write_text(json.dumps(z,indent=2,allow_nan=False)+'\n')
    temp.replace(a.state)
    print(json.dumps(r,indent=2,allow_nan=False))
