"""Opt-in continuation with independent Gauss scalar nodal occupations."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_scalar_gauss2_live_fixture import ScalarGauss2
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load
from early_pole_resolved_parent_reference import S


def run(g,div,s,cut,floor,target,n,restart=None):
    d=ScalarGauss2(load(g)[0],SPECS[g],div,s,cut,gauss_floor=floor)
    key=dict(grid=g,charge_div=div,scalar_div=s,cut=cut,gauss_floor=floor,reference_x=X,
             nodes=[d.ns,d.nf,d.np],scalar_cell_count=d.old_ns,scalar_nodal_rule='gauss-volume-v1')
    x0=X
    init=None
    if restart is not None:
        p=json.loads(restart.read_text())
        if any(p.get(k)!=v for k,v in key.items()):
            raise ArithmeticError('Restart changes the numerical nodal chart')
        x0=p['x_end']
        init=np.asarray(p['nodal_occupation'])
        if init.shape!=(d.n,) or not np.all(np.isfinite(init)) or min(init)<0 or max(init[d.ns:])>1:
            raise ArithmeticError('Invalid saved Gauss nodal state')
    if not X<=x0<target or n<8 or n%8:
        raise ValueError('Invalid clock or mesh count')
    A,_=drive(d,target,n,x_start=x0,initial=init,stepper=step)
    B,end=drive(d,target,2*n,x_start=x0,initial=init,stepper=step)
    drift=abs(A[-1]['scalar_yield']/B[-1]['scalar_yield']-1)
    if drift>1e-3 or max(abs(t['continuity_residual']) for t in A+B)>3e-10:
        raise ArithmeticError('Gauss history has not passed its time or continuity gate')
    pos=np.r_[0,np.cumsum(d.scalar_repeats)[:-1]]
    per_cell=np.add.reduceat(end[:d.ns]*d.w[:d.ns],pos)/S
    if abs(fsum(per_cell)/B[-1]['scalar_yield']-1)>4e-14:
        raise ArithmeticError('Scalar Gauss cell counts do not close')
    card=dict(key,x_end=target,steps=2*n,nodal_occupation=end.tolist())
    report=dict(key,starting_x=x0,x_end=target,coarse_steps=n,relative_time_drift=drift,
                checkpoint_trajectories=[A,B],scalar_cell_yields=per_cell.tolist(),
                scalar_pole_mass=d.state(target,end)[0],
                scope='Alternate scalar nodal quadrature; do not interpret nodal occupations as constant within cells.')
    return report,card


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=SPECS)
    p.add_argument('--charge-div',type=int,default=2)
    p.add_argument('--scalar-div',type=int,default=4)
    p.add_argument('--cut',type=float,default=1e-5)
    p.add_argument('--gauss-floor',type=float)
    p.add_argument('--x-end',type=float,default=.0075)
    p.add_argument('--steps',type=int,default=2048)
    p.add_argument('--state',type=Path)
    p.add_argument('--restart',type=Path)
    a=p.parse_args()
    r,s=run(a.grid,a.charge_div,a.scalar_div,a.cut,a.gauss_floor,a.x_end,a.steps,a.restart)
    if a.state:
        tmp=a.state.with_suffix(a.state.suffix+'.tmp')
        tmp.write_text(json.dumps(s,indent=2,allow_nan=False)+'\n')
        tmp.replace(a.state)
    print(json.dumps(r,indent=2,allow_nan=False))
