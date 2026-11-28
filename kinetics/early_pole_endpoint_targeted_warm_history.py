"""Separate evidence path for the predeclared targeted positive warm chart."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_chain_receipt import verify_chart,CLOCK,EPOCHS,read
from early_pole_endpoint_higher_warm_gauss_history import run_chart
from early_pole_endpoint_stencil_remap_pilot import X,load
from early_pole_endpoint_warm_targeted_source_refinement import PARENTS,Targeted


def build():
    cfg=load('fine')[0]
    d=Targeted(cfg)
    picked=d.ds.ed[:-1][d.repeats==8].tolist()
    if picked!=list(PARENTS) or np.count_nonzero(d.repeats==4)!=20 or (d.ns,d.nf,d.np)!=(219,448,488):
        raise ArithmeticError('Targeted test chart changed')
    if np.min(d.ds.w)<=0 or not np.all(np.isfinite([d.ds.p,d.ds.w])):
        raise ArithmeticError('Nonpositive or nonfinite scalar volume')
    n=np.r_[0,np.cumsum(d.repeats)]
    q=np.repeat(d.ds.ed[:-1],d.repeats)
    r=np.repeat(d.ds.ed[1:],d.repeats)
    if np.any(d.ds.p<=q) or np.any(d.ds.p>=r) or len(d.initial)!=d.n or n[-1]!=d.ns:
        raise ArithmeticError('Lifted initial state or nodes do not align')
    meta=dict(grid='fine', charge_div=2, scalar_div=4, cut=1e-5, gauss_floor=.01,
              reference_x=X,nodes=[d.ns,d.nf,d.np],scalar_cell_count=d.old_count,
              scalar_nodal_rule='gauss-volume-v3-targeted',
              targeted_rule=dict(baseline_gauss_order=4,target_q_lower_edges=list(PARENTS),
                                 split_fraction=.5,gauss_order_per_half=4))
    return d,meta


def run(end,steps,restart):
    d,m=build()
    if restart is not None:
        u,_=read(restart)
        if any(u.get(k)!=v for k,v in m.items()):
            raise ArithmeticError('Mismatched targeted warm restart')
        ix=EPOCHS.index(u['x_end'])
        if ix==len(EPOCHS)-1 or u['steps']!=2*CLOCK[ix] or (end,steps)!=(EPOCHS[ix+1],CLOCK[ix+1]):
            raise ArithmeticError('Targeted restart does not select the unique next clock segment')
    return run_chart(d,m,end,steps,restart,scope='Independent targeted positive warm history; fixed cold and charge charts, not a physical abundance result.')


def audit(reports,states,*,prefix=False):
    d,m=build()
    n=len(reports)
    if n!=len(states) or not 1<=n<=len(EPOCHS) or not prefix and n!=len(EPOCHS):
        raise ArithmeticError('Must supply paired ordered files for the full chain or an explicit prefix')
    out=[]
    prior=None
    for i,(p,q) in enumerate(zip(reports,states,strict=True)):
        r,rh=read(p); s,sh=read(q)
        row=verify_chart(r,s,d,(X,*EPOCHS[:-1])[i],EPOCHS[i],CLOCK[i],m,d.repeats)
        if any('actual_previous_raw_sha256' not in z or z['actual_previous_raw_sha256']!=prior for z in (r,s)):
            raise ArithmeticError('Targeted chain did not use the preceding uncompressed fine bytes')
        out.append(dict(**row,verified_previous_raw_sha256=prior,source_report_sha256=rh,uncompressed_state_sha256=sh))
        prior=sh
    return dict(complete=n==len(EPOCHS),**m,checked_epoch_pairs=out,
                scope='Check actual paired trajectory archives and restart links, not a replay of the ODE.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='mode',required=True)
    a=sub.add_parser('advance')
    a.add_argument('--x-end',type=float,default=.0075)
    a.add_argument('--steps',type=int,default=2048)
    a.add_argument('--restart',type=Path)
    a.add_argument('--state',type=Path,required=True)
    c=sub.add_parser('audit')
    c.add_argument('--reports',type=Path,nargs='+',required=True)
    c.add_argument('--states',type=Path,nargs='+',required=True)
    c.add_argument('--prefix',action='store_true')
    t=p.parse_args()
    if t.mode=='advance':
        r,z=run(t.x_end,t.steps,t.restart)
        tmp=t.state.with_suffix('.part')
        tmp.write_text(json.dumps(z,indent=2,allow_nan=False)+'\n')
        tmp.replace(t.state)
    else:
        r=audit(t.reports,t.states,prefix=t.prefix)
    print(json.dumps(r,indent=2,allow_nan=False))
