"""Impose a scalar 4-to-2 count change at fixed scale, without assigning a rate."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_decay_mode_clock_reference import Q, intervals
from early_pole_massive_daughter_energy_stencil_reference import MassiveFermion


def stress(d,x,f):
    n=d.d.w*f
    mb=1e-9*x
    m=brentq(lambda u:u*u-mb*mb-.5*d.d.lam*(n@(1/np.hypot(u,d.d.p))),
             mb,mb+1.01*(.5*d.d.lam*n.sum())**(1/3),xtol=2e-19,rtol=1e-15)
    assert abs(m/d.d.gap(x,f)-1)<2e-11
    E=np.hypot(m,d.d.p)
    U=(m*m-mb*mb)**2/(2*d.d.lam)
    rho=E@n-U
    pressure=(d.d.p**2/(3*E))@n+U
    trace=mb*mb*(n@(1/E))
    assert abs((rho-3*pressure)-trace)<3e-14*rho
    return m,E,rho,pressure,trace


def width(x,b,m,p,pe,fe,P,B):
    E,l,h,span,_=intervals(x,m,b,p)
    assert span>0. and abs((h-l)/span-1)<1e-10
    pp=np.hypot(x,pe)
    df=E+np.hypot(b,fe)
    cut=np.unique(np.r_[l,h, pp[(pp>l)&(pp<h)],df[(df>l)&(df<h)]])
    w=np.diff(cut); c=.5*(cut[1:]+cut[:-1])
    i=np.searchsorted(pp,c)-1; j=np.searchsorted(df,c)-1
    assert min(i)>=0 and min(j)>=0 and max(i)<len(P) and max(j)<len(B)
    u,v=np.asarray(P)[i],np.asarray(B)[j]
    return dict(daughter_energy_center_over_T=(l+h)/2-E,original_half_span_over_T=span/2,
                average_parent_occupation=float(w@u/span),
                average_daughter_occupation=float(w@v/span),
                average_net_stimulus=float(w@(u-v)/span),
                daughter_shells=len(np.unique(j)))


def impose(d,s,loss,steps):
    f=np.asarray(s['endpoint_scalar_occupations']).copy()
    f=np.maximum(f,0.)
    N=d.d.w*f
    n0=N.sum()
    q=d.d.p
    x=s['x_end']
    m,E,rho,pres,tr=stress(d,x,f)
    initial=lambda mu:width(x,d.ratio*x,mu,Q,d.pe,d.d.fed,
                            s['endpoint_parent_occupations'],s['endpoint_fermion_occupations'])
    peak=initial(m)
    cold=d.d.ed[1:]<=1e-5
    counts=loss*n0/(2*steps)
    maxerr=0.
    for _ in range(steps):
        mm,ee,rr,_,_=stress(d,x,N/d.d.w)
        remove=4*counts*N[cold]/N[cold].sum()
        def event(energies):
            target=(remove@energies[cold])/(2*counts)
            j=np.searchsorted(energies,target)
            assert 0<j<len(q)
            a=(target-energies[j-1])/(energies[j]-energies[j-1])
            dn=np.zeros_like(N)
            dn[cold]-=remove
            dn[j-1]+=2*counts*(1-a)
            dn[j]+=2*counts*a
            return dn
        euler=event(ee)
        mm,ee,_,_,_=stress(d,x,(N+.5*euler)/d.d.w)
        dn=event(ee)
        err=abs(dn@ee)/(4*counts*mm)
        maxerr=max(maxerr,err)
        N+=dn
        assert np.min(N)>=0.
    m2,E2,r2,p2,t2=stress(d,x,N/d.d.w)
    end=initial(m2)
    assert abs(N.sum()/n0-(1-loss))<1e-12
    return dict(loss=loss,steps=steps,initial_pole=m,final_pole=m2,
                dlog_pole_over_dlog_scalar_count=np.log(m2/m)/np.log(1-loss),
                energy_relative_error=(r2/rho-1),max_event_energy_error=maxerr,
                initial_w=pres/rho,final_w=p2/r2,initial_trace=tr,final_trace=t2,
                mode_q=Q,window_initial=peak,window_after=end,
                window_center_shift_in_initial_half_spans=
                    (end['daughter_energy_center_over_T']-peak['daughter_energy_center_over_T'])/peak['original_half_span_over_T'])


def one(s):
    shape={(67,112,122):(65,90,72),(97,156,176):(95,125,108)}
    d=MassiveFermion(*shape[tuple(s['cells'])],s['mass_ratio'],0.,True)
    rows=[]
    for loss in (.0001,.01,.05):
        for steps in (16,128):
            rows.append(impose(d,s,loss,steps))
    assert abs(rows[1]['dlog_pole_over_dlog_scalar_count']-.5)<.004
    assert all(abs(row['energy_relative_error'])<3e-8 for row in rows)
    for a,b in zip(rows[::2],rows[1::2],strict=True):
        assert abs(b['energy_relative_error'])<abs(a['energy_relative_error'])/5+1e-13
    return dict(x=s['x_end'],mass_ratio=s['mass_ratio'],results=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    a=p.parse_args()
    out={}
    for name in a.spectra:
        f=Path(name)
        raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        out[f.name]=one(json.loads(raw))
    print(json.dumps(out,indent=2,allow_nan=False))
