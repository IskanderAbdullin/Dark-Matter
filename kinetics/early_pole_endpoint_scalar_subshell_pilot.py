"""Opt-in scalar-only resolution with fixed charge chart and unchanged clock."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load
from early_pole_resolved_parent_reference import S


class ScalarSubshellChart(CellChart):
    """Keep the native scalar COM convention and split only shells inside the cold cut."""
    def __init__(self,bundle,shape,charge_div,scalar_div,cut,*,warm_floor=None,warm_div=1):
        super().__init__(bundle,shape,charge_div,reference_x=X)
        edges=self.ds.ed
        split=np.where(edges[1:]<=cut,scalar_div,1)
        if scalar_div<1 or not np.any(edges[1:]<=cut):
            raise ValueError('No eligible cold scalar shells or invalid scalar refinement')
        if warm_floor is not None:
            warm=edges[:-1]>=warm_floor
            if not np.isfinite(warm_floor) or warm_floor<cut or warm_div<1 or not np.any(warm):
                raise ValueError('Invalid or overlapping warm scalar range')
            split=np.where(warm,warm_div,split)
        elif warm_div!=1:
            raise ValueError('Warm subdivisions require an explicit warm floor')
        # Endpoints remain exactly those of the archived chart.
        children=[np.geomspace(lo,hi,n+1)[:-1] for lo,hi,n in zip(edges[:-1],edges[1:],split,strict=True)]
        self.ds.ed=np.r_[np.concatenate(children),edges[-1]]
        if not np.array_equal(self.ds.ed[np.r_[0,np.cumsum(split)]],edges):
            raise ArithmeticError('A native shell edge moved')
        self.ds.p=np.sqrt(self.ds.ed[:-1]*self.ds.ed[1:])
        self.ds.w=np.diff(self.ds.ed**3)/(6*np.pi**2)
        self.ds.n=len(self.ds.w)
        parent_cell=np.repeat(np.arange(self.ns),split)
        old=self.w[:self.ns]*self.initial[:self.ns]
        new=np.repeat(self.initial[:self.ns],split)
        inherited=np.bincount(parent_cell,weights=self.ds.w*new,minlength=len(old))
        if np.max(np.abs(inherited-old)/np.maximum(old,1e-70))>3e-14:
            raise ArithmeticError('Lost native scalar shell number in control remap')
        self.initial=np.r_[new,self.initial[self.ns:]]
        self.ns=self.ds.n
        self.iP=self.ns+self.nf
        self.n=self.iP+self.np
        self.w=np.r_[self.ds.w,self.wd,self.wp]
        self.seed_origin[0]=self.ds.ed
        self.native_scalar_labels=parent_cell


def run(grid,charge_div,scalar_div,cut,xend,n,restart=None,fast_capacity=False,warm_floor=None,warm_div=1):
    a=load(grid)[0]
    chart=ScalarSubshellChart(a,SPECS[grid],charge_div,scalar_div,cut,warm_floor=warm_floor,warm_div=warm_div)
    key=dict(grid=grid,charge_div=charge_div,scalar_div=scalar_div,cut=cut,
             counts=[chart.ns,chart.nf,chart.np],reference_x=X)
    if warm_floor is not None:
        key.update(warm_floor=warm_floor,warm_div=warm_div)
    x0=X
    y=None
    if restart:
        prior=json.loads(restart.read_text())
        if any(prior.get(k)!=v for k,v in key.items()) or prior.get('warm_floor')!=warm_floor:
            raise ArithmeticError('Changed scalar/charge chart in restart')
        x0=prior['x_end']
        y=np.array(prior['endpoint_occupation'])
        if y.shape!=(chart.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[chart.ns:])>1:
            raise ArithmeticError('Invalid restart state')
    if not X<=x0<xend or n<8 or n%8:
        raise ValueError('Invalid clock or sample count')
    kw={}
    if fast_capacity:
        from early_pole_endpoint_capacity_fast_step import step
        kw['stepper']=step
    low,_=drive(chart,xend,n,x_start=x0,initial=y,**kw)
    high,z=drive(chart,xend,2*n,x_start=x0,initial=y,**kw)
    check=abs(low[-1]['scalar_yield']/high[-1]['scalar_yield']-1)
    if check>1e-3 or max(abs(r['continuity_residual']) for r in low+high)>3e-10:
        raise ArithmeticError(('Failed temporal refinement or continuity ledger',check))
    state=dict(key,x_end=xend,steps=2*n,endpoint_occupation=z.tolist())
    native_num=np.bincount(chart.native_scalar_labels,weights=chart.ds.w*z[:chart.ns])
    report=dict(key,starting_x=x0,x_end=xend,coarse_steps=n,relative_time_drift=check,
                checkpoint_trajectories=[low,high],native_scalar_yield=(native_num/S).tolist(),
                scalar_pole_mass=chart.state(xend,z)[0],
                scope='Scalar subshell isolation at fixed charge chart, not a physical spectrum.')
    if abs(fsum(report['native_scalar_yield'])/high[-1]['scalar_yield']-1)>4e-14:
        raise ArithmeticError('Native aggregation changed scalar count')
    return report,state


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=SPECS)
    p.add_argument('--charge-div',type=int,default=2)
    p.add_argument('--scalar-div',type=int,default=2)
    p.add_argument('--cut',type=float,default=1e-5)
    p.add_argument('--x-end',type=float,default=.0075)
    p.add_argument('--steps',type=int,default=1024)
    p.add_argument('--state',type=Path)
    p.add_argument('--restart',type=Path)
    p.add_argument('--fast-capacity',action='store_true')
    p.add_argument('--warm-floor',type=float)
    p.add_argument('--warm-div',type=int,default=1)
    o=p.parse_args()
    r,s=run(o.grid,o.charge_div,o.scalar_div,o.cut,o.x_end,o.steps,o.restart,o.fast_capacity,o.warm_floor,o.warm_div)
    if o.state:
        t=o.state.with_suffix('.tmp')
        t.write_text(json.dumps(s,indent=2,allow_nan=False)+'\n')
        t.replace(o.state)
    print(json.dumps(r,indent=2,allow_nan=False))
