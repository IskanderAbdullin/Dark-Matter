"""Step-doubling while accepting only bounded projected reaction updates."""
import argparse
import json
from math import fsum
import sys

import numpy as np

from early_pole_bounded_face_extents_pilot import step
from early_pole_double_charge_stencil_pilot import chart
from early_pole_edge_padded_charged_projection_pilot import charge_mesh
from early_pole_endpoint_stencil_remap_pilot import LABELS,SPECS,archives,load,parts,remap,MassiveFermion
from early_pole_live_bounded_face_extents_pilot import energy_difference,hold_trace
from early_pole_live_parent_energy_reference import state_moments


def once(d,x0,state,t,h):
    x=x0*np.exp(t+h/2)
    _,ra,rb,m=state_moments(d,x,state)
    H=np.sqrt(d.radiation/(d.radiation+ra+rb))
    C,_=chart(d,x,m,4)
    half,l=step(C,d,state,H*h/2)
    _,aa,bb,mm=state_moments(d,x,half)
    U=np.sqrt(d.radiation/(d.radiation+aa+bb))
    B,_=chart(d,x,mm,4)
    end,k=step(B,d,state,U*h,rate_state=half)
    if min(end)<0 or max(end[d.ns:])>1:
        raise ArithmeticError('An accepted candidate would violate the occupation caps')
    midpoint=t+h/2
    trace=hold_trace(d,x0,t,midpoint,state)+hold_trace(d,x0,midpoint,t+h,end)
    jump=energy_difference(d,x,state,x,end)
    return end,max(l,k),trace,jump


def drive(d,x0,y,span,tol):
    t=0.
    h=min(2e-7,span)
    trigger=0
    rejected=0
    x=y.copy()
    ts,js,dt=[],[],[]
    maxf=max(x[d.ns:])
    lim=0.
    worst=0.
    while t<span and len(dt)<20000:
        h=min(h,span-t)
        if h<1e-13 and span-t>1e-13:
            raise ArithmeticError('Unresolved scales before reaching the planned adaptive end')
        a,*_=once(d,x0,x,t,h)
        m,c,dtrace,j=once(d,x0,x,t,h/2)
        b,k,etrace,q=once(d,x0,m,t+h/2,h/2)
        w=d.w
        diff=np.split(w*abs(b-a),(d.ns,d.iP))
        move=np.split(w*abs(b-x),(d.ns,d.iP))
        # Weighted L1 error detects redistribution despite cancellation
        # of the three species' total numbers. A separate pointwise cap
        # protects appreciably changing charged occupations.
        l1=max(fsum(r)/(tol*max(fsum(f),h*1e-9)) for r,f in zip(diff,move,strict=True))
        charged=np.max(abs(b[d.ns:]-a[d.ns:])/(2e-4+4*tol*abs(b[d.ns:]-x[d.ns:])))
        er=max(l1,charged)
        if er<=1:
            x=b
            t+=h
            dt.append(h)
            ts.extend([dtrace,etrace])
            js.extend([j,q])
            lim=max(lim,c,k)
            maxf=max(maxf,max(x[d.ns:]))
            worst=max(worst,er)
        else:
            rejected+=1
        h*=min(2.,max(.2,.8/np.sqrt(max(er,1e-16))))
        trigger+=1
        if trigger%128==0:
            print(f'accepted={len(dt)} rejected={rejected} t={t:.12g} h={h:.6g}',file=sys.stderr,flush=True)
    if t<span:
        raise ArithmeticError('Incomplete adaptive interval')
    w=d.w*(x-y)
    dn=[float(fsum(v)) for v in np.split(w,(d.ns,d.iP))]
    e=energy_difference(d,x0,y,x0*np.exp(t),x)
    residual=e-fsum(ts)
    count=max(abs(dn[0]-dn[1]),abs(dn[0]+dn[2]))
    if count>3e-12 or abs(residual-fsum(js))>2e-14:
        raise ArithmeticError(('Independent adaptive budgets do not close',count,residual,fsum(js)))
    return dict(span=t,tol=tol,accepted=len(dt),rejected=rejected,
                min_max_accepted_step=[min(dt),max(dt)],max_error_estimate=worst,
                paired_count_changes=dn,number_defect=count,maximum_charged_occupation=float(maxf),
                maximum_limited_extent_fraction=lim,continuity_residual=residual,
                integrated_collision_energy_jumps=fsum(js))


def run(g,span,tolerances):
    a=load(g)[0]
    x=a['x_end']
    d=MassiveFermion(*SPECS[g],a['mass_ratio'],a['elastic_rate'],True)
    charge_mesh(d,1/4096,1/4096)
    y=np.concatenate([remap(o,a['endpoint_'+k+'_occupations'],dest)[0]
                      for o,k,dest in zip(archives(a),LABELS,parts(d),strict=True)])
    return dict(grid=g,scope='Archive-started bounded adaptive pilot',cases=[drive(d,x,y,span,tol) for tol in tolerances])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=SPECS)
    p.add_argument('--span',type=float,default=5e-5)
    p.add_argument('--tol',nargs='+',type=float,default=[2e-3,5e-4])
    a=p.parse_args()
    if a.span<=0 or min(a.tol)<=0:
        raise ValueError('Positive span and error thresholds required')
    print(json.dumps(run(a.grid,a.span,a.tol),indent=2,allow_nan=False))
