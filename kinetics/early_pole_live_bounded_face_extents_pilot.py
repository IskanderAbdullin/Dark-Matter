"""Resolve an archive-started FRW increment with bounded reaction extents."""
import argparse
import json
from math import fsum

import numpy as np

from early_pole_bounded_face_extents_pilot import step
from early_pole_double_charge_stencil_pilot import chart
from early_pole_edge_padded_charged_projection_pilot import charge_mesh
from early_pole_endpoint_stencil_remap_pilot import LABELS,SPECS,X,archives,load,parts,remap,MassiveFermion
from early_pole_live_parent_energy_reference import state_moments


def hold_trace(d,x0,lo,hi,f):
    m=(lo+hi)/2
    w=(hi-lo)/2
    n=w/np.sqrt(3)
    t=np.log(x0/X)
    return w*(d.trace(t+m-n,f)+d.trace(t+m+n,f))



def energy_difference(d,x0,y0,x1,y1):
    # Difference first: subtracting two O(1) total energies masks the
    # O(1e-10) collision signal in these short replays.
    m0=d.d.gap(x0,y0[:d.ns])
    m1=d.d.gap(x1,y1[:d.ns])
    ea=np.hypot(d.d.p,m0)
    eb=np.hypot(d.d.p,m1)
    dEs=(m1-m0)*(m1+m0)/(ea+eb)
    A=np.hypot(d.daughter_q,d.daughter_mass(x0))
    B=np.hypot(d.daughter_q,d.daughter_mass(x1))
    u=d.daughter_mass(x0)
    v=d.daughter_mass(x1)
    dEd=np.sum(d.daughter_w*((v-u)*(v+u)/(A+B)),axis=1)/d.d.wf
    P=np.hypot(d.pq,x0)
    Q=np.hypot(d.pq,x1)
    dEp=np.sum(d.pqw*((x1-x0)*(x1+x0)/(P+Q)),axis=1)/d.pw
    E0=np.r_[ea,d.daughter_energies(x0),d.parent_energies(x0)]
    dE=np.r_[dEs,dEd,dEp]
    b0=1e-9*x0
    b1=1e-9*x1
    dx=(m1-m0)*(m1+m0)-(b1-b0)*(b1+b0)
    dU=dx*(m0*m0+m1*m1-b0*b0-b1*b1)/(2*d.d.lam)
    return fsum(d.w*((y1-y0)*E0+y1*dE))-dU

def trial(d,x0,z,span,n,midpoint):
    f=z.copy()
    t=np.linspace(0,span,n+1)
    trace=[]
    jumps=[]
    lim=[]
    maxf=0.
    for lo,hi in zip(t[:-1],t[1:],strict=True):
        s=(lo+hi)/2
        x=x0*np.exp(s)
        er,a,b,m=state_moments(d,x,f)
        H=np.sqrt(d.radiation/(d.radiation+a+b))
        ch,_=chart(d,x,m,4)
        if midpoint:
            predictor,l0=step(ch,d,f,H*(hi-lo)/2)
            _,aa,bb,mm=state_moments(d,x,predictor)
            Hm=np.sqrt(d.radiation/(d.radiation+aa+bb))
            pair,_=chart(d,x,mm,4)
            g,limiter=step(pair,d,f,Hm*(hi-lo),rate_state=predictor)
            limiter=max(limiter,l0)
        else:
            g,limiter=step(ch,d,f,H*(hi-lo))
        if min(g)<0 or np.max(g[d.ns:])>1:
            raise ArithmeticError('A live substep left the occupation bounds')
        # Integrate the trace of each collisionless half at independent
        # Gaussian nodes. Also record the exact energy jump at each collision.
        trace.append(hold_trace(d,x0,lo,s,f)+hold_trace(d,x0,s,hi,g))
        jumps.append(energy_difference(d,x,f,x,g))
        f=g
        maxf=max(maxf,np.max(f[d.ns:]))
        lim.append(limiter)
    E=energy_difference(d,x0,z,x0*np.exp(span),f)
    excess=E-fsum(trace)
    h=(f-z)*d.w
    N=[fsum(a) for a in np.split(h,(d.ns,d.iP))]
    delta=max(abs(N[0]-N[1]),abs(N[0]+N[2]))
    if delta>6e-14 or abs(excess-fsum(jumps))>5e-16:
        raise ArithmeticError(('Live number or independently integrated trace control',delta,excess,fsum(jumps)))
    return dict(steps=n,paired_delta_N=N,number_error=delta,
                continuity_residual=excess,sum_fixed_clock_collision_energy_jumps=fsum(jumps),
                delta_total_energy=E,independent_collisionless_trace=fsum(trace),
                maximum_charged_occupation=maxf,maximum_limited_extent_fraction=max(lim),
                scalar_pole_at_end=state_moments(d,x0*np.exp(span),f)[-1])


def run(g,span,midpoint=False,base=32):
    a=load(g)[0]
    x0=a['x_end']
    d=MassiveFermion(*SPECS[g],a['mass_ratio'],a['elastic_rate'],True)
    charge_mesh(d,1/4096,1/4096)
    values=[a['endpoint_'+s+'_occupations'] for s in LABELS]
    z=np.concatenate([remap(o,f,e)[0] for o,f,e in zip(archives(a),values,parts(d),strict=True)])
    cases=[trial(d,x0,z,span,n,midpoint) for n in (base,2*base,4*base)]
    a,b,c=[r['paired_delta_N'][0] for r in cases]
    diff=float(abs(c-b)/max(abs(c),1e-50))
    if diff>2e-3 or abs(c-b)>=abs(b-a):
        raise ArithmeticError(('Live time refinement not resolved for the declared interval',g,span,a,b,c))
    return dict(grid=g,span_in_log_scale_factor=span,time_scheme='midpoint' if midpoint else 'euler',scope='Local archive-started pilot, not a full time history.',
                last_relative_scalar_number_change=diff,refinement_ratio=float((b-a)/(c-b)),data=cases)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=SPECS)
    p.add_argument('--span',type=float,default=1e-6)
    p.add_argument('--midpoint',action='store_true')
    p.add_argument('--base',type=int,default=32)
    a=p.parse_args()
    print(json.dumps(run(a.grid,a.span,a.midpoint,a.base),indent=2,allow_nan=False))
