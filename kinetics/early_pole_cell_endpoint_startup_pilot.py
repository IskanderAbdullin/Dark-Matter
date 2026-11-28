"""Run opt-in endpoint cells from an actual empty-daughter startup."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import expit

from early_pole_bounded_face_extents_pilot import step
from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load
from early_pole_live_bounded_face_extents_pilot import hold_trace
from early_pole_resolved_parent_reference import S


def checkpoint_signature(d,y,t,x_start,xend,n,stepper,provenance):
    if (not isinstance(provenance,str) or len(provenance)!=64
            or set(provenance)-set('0123456789abcdef') or n<8 or not 0.<x_start<xend):
        raise ValueError('Resumable drive requires locked provenance and a valid positive clock')
    a=hashlib.sha256()
    a.update(json.dumps([provenance,stepper.__module__,stepper.__name__,x_start,xend,n,d.ns,d.iP,d.n]).encode())
    for z in (d.ds.ed,d.ds.p,d.ds.w,d.wp,d.wd,d.qp,d.qd,y,t):
        b=np.asarray(z,dtype='<f8')
        a.update(str(b.shape).encode()+b.tobytes())
    return a.hexdigest()


def drive_checkpoint_digest(packet):
    payload={k:v for k,v in packet.items() if k!='payload_sha256'}
    return hashlib.sha256(json.dumps(payload,allow_nan=False,separators=(',',':'),sort_keys=True).encode()).hexdigest()


def write_drive_checkpoint(path,packet):
    target=Path(path)
    tmp=target.with_suffix(target.suffix+'.tmp')
    p=dict(packet,payload_sha256=drive_checkpoint_digest(packet))
    tmp.write_text(json.dumps(p,allow_nan=False,separators=(',',':'))+'\n')
    tmp.replace(target)


def drive(d,xend,n,*,x_start=X,initial=None,stepper=step,checkpoint=None,checkpoint_stride=0,
          restart=None,provenance=None):
    folds=len(d.wp)//(2*d.d0.np)
    parent=np.repeat(expit(-d.d0.parent_energies(X)),2*folds)
    y=np.r_[np.zeros(d.iP),parent] if initial is None else initial.copy()
    state=y.copy()
    t=np.linspace(0.,np.log(xend/x_start),n+1)
    trace=[]; caps=[]; checks=[]
    cursor=0
    key=None
    if checkpoint is not None or restart is not None:
        key=checkpoint_signature(d,y,t,x_start,xend,n,stepper,provenance)
        if checkpoint is not None and (not callable(checkpoint) or type(checkpoint_stride) is not int or checkpoint_stride<=0):
            raise ValueError('Checkpoint requires a callback and positive integral stride')
        if restart is not None:
            if restart.get('payload_sha256')!=drive_checkpoint_digest(restart):
                raise ArithmeticError('Corrupt step-checkpoint payload')
            if restart.get('signature')!=key:
                raise ArithmeticError('Step checkpoint does not bind this exact clock, initial state and geometry')
            cursor=restart.get('completed')
            trace=list(restart['trace'])
            caps=list(restart['caps'])
            checks=list(restart['checks'])
            state=np.asarray(restart['occupation'],dtype=float)
            if (type(cursor) is not int or not 1<=cursor<=n or len(trace)!=cursor or len(caps)!=cursor
                    or len(checks)!=cursor//(n//8)+int(cursor==n and cursor%(n//8)!=0)
                    or restart.get('target_steps')!=n or restart.get('x_reached')!=float(x_start*np.exp(t[cursor]))
                    or state.shape!=(d.n,) or not np.all(np.isfinite(np.r_[state,trace,caps]))
                    or min(state)<0. or max(state[d.ns:])>1.):
                raise ArithmeticError('Invalid checkpoint state or unresolved continuity ledger')
    for index,(a,b) in enumerate(zip(t[cursor:-1],t[cursor+1:],strict=True),cursor):
        mid=(a+b)/2
        x=x_start*np.exp(mid)
        M,_,_,H=d.state(x,state)
        c,_=d.chart(x,M,3)
        predictor,cap1=stepper(c,d,state,H*(b-a)/2)
        M,_,_,H=d.state(x,predictor)
        c,_=d.chart(x,M,3)
        now,cap2=stepper(c,d,state,H*(b-a),rate_state=predictor)
        trace.append(hold_trace(d,x_start,a,mid,state)+hold_trace(d,x_start,mid,b,now))
        state=now
        caps.append(max(cap1,cap2))
        if (index+1)%(n//8)==0 or index+1==n:
            s,f,p=[fsum(w) for w in np.split(d.w*(state-y),(d.ns,d.iP))]
            e=d.energy_difference(x_start,y,x_start*np.exp(b),state)-fsum(trace)
            if min(state)<0 or max(state[d.ns:])>1 or max(abs(s-f),abs(s+p))>3e-13:
                raise ArithmeticError('Startup Pauli or paired number gate')
            checks.append(dict(x=float(x_start*np.exp(b)),scalar_yield=fsum(d.w[:d.ns]*state[:d.ns])/S,pair_number_error=max(abs(s-f),abs(s+p)),
                               continuity_residual=e,max_charged=float(max(state[d.ns:])),
                               max_limited_extent=float(max(caps))))
        if checkpoint is not None and ((index+1)%checkpoint_stride==0 or index+1==n):
            record=dict(signature=key,completed=index+1,target_steps=n,x_reached=float(x_start*np.exp(b)),
                        occupation=state.tolist(),trace=trace.copy(),
                        caps=caps.copy(),checks=[dict(z) for z in checks])
            checkpoint(dict(record,payload_sha256=drive_checkpoint_digest(record)))
    return checks,state


def run(g,xend,div,n,state_path=None,one=False,restart=None,fast_capacity=False):
    a=load(g)[0]
    d=CellChart(a,SPECS[g],div,reference_x=X)
    x0=X
    initial=None
    if restart is not None:
        saved=json.loads(Path(restart).read_text())
        if saved['grid']!=g or saved['subcells_per_archived_cell']!=div or saved['reference_x']!=X or saved['counts']!=[d.ns,d.nf,d.np]:
            raise ArithmeticError('Restart does not belong to this endpoint chart')
        x0=saved['x_end']
        initial=np.array(saved['endpoint_occupation'])
        if not X<x0<xend or initial.shape!=(d.n,) or not np.all(np.isfinite(initial)) or min(initial)<0 or max(initial[d.ns:])>1:
            raise ArithmeticError('Invalid restart interval or occupations')
    stepper=step
    if fast_capacity:
        from early_pole_endpoint_capacity_fast_step import step as stepper
    arm=lambda k: drive(d,xend,k,x_start=x0,initial=initial,stepper=stepper)
    if one:
        hard,hard_state=arm(n)
        result=[hard]
        n_high=n
    else:
        coarse,_=arm(n)
        hard,hard_state=arm(2*n)
        result=[coarse,hard]
        n_high=2*n
    if state_path is not None:
        card=dict(grid=g,subcells_per_archived_cell=div,counts=[d.ns,d.nf,d.np],x_end=xend,
                  steps=n_high,reference_x=X,statistics='constant native per-bin parent seeds, empty daughters',
                  endpoint_occupation=hard_state.tolist())
        f=Path(state_path)
        tmp=f.with_suffix(f.suffix+'.tmp')
        tmp.write_text(json.dumps(card,indent=2,allow_nan=False)+'\n')
        tmp.replace(f)

    drift=abs(result[0][-1]['scalar_yield']/result[-1][-1]['scalar_yield']-1) if not one else None
    old=[c['scalar_yield'] for c in a['checkpoints'] if c['x']==xend]
    output=dict(grid=g,subcells_per_archived_cell=div,counts=[d.ns,d.nf,d.np],x_end=xend,
                time_refinement=n,endpoint_yield_relative_drift=drift,
                original_archive_yield=old[0] if old else None,
                checkpoint_trajectories=result,
                scope='Opt-in startup check; cannot stand in for the full archived production history')
    if restart is not None:
        output.update(starting_x=x0,starting_scalar_yield=fsum(d.w[:d.ns]*initial[:d.ns])/S)
    return output


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=SPECS)
    p.add_argument('--x-end',type=float,default=.0075)
    p.add_argument('--div',type=int,default=1)
    p.add_argument('--steps',type=int,default=256)
    p.add_argument('--state-path',type=Path)
    p.add_argument('--one',action='store_true')
    p.add_argument('--fast-capacity',action='store_true',help='Opt in to the proved no-limiting fast path')
    p.add_argument('--restart',type=Path)
    a=p.parse_args()
    if not X<a.x_end or a.steps<8 or a.div<1:
        raise ValueError('Invalid startup interval or resolution')
    print(json.dumps(run(a.grid,a.x_end,a.div,a.steps,a.state_path,a.one,a.restart,a.fast_capacity),indent=2,allow_nan=False))
