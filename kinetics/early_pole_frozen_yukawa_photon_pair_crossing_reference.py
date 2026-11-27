"""Crossed real-photon channel: photon and scalar versus two fermion flavors."""
import argparse
import gzip
import json

import numpy as np
from scipy.special import roots_legendre

from early_pole_hard_yukawa_photoproduction_reference import amplitudes,invariant_form
from early_pole_frozen_yukawa_photoproduction_source_reference import edges_between
from early_pole_native_aperture_witness_reference import GRID,occup


def angular(s,c,M,m,ms):
    root=np.sqrt(s)
    pr=np.sqrt((s-(M+m)**2)*(s-(M-m)**2))/(2*root)
    Es=np.hypot(pr,M)
    Eb=np.hypot(pr,m)
    k=(s-ms*ms)/(2*root)
    a=-k*(M*M/(Es+pr)+pr*(1+c))
    b=k*(m*m/(Eb+pr)+pr*(1-c))
    d=-(s-M*M-m*m)/2
    return -invariant_form(a,b,d,M,m)


def trace_control(M,m,ms):
    c=np.r_[-.998,-.8,0.,.8,.998]
    errors=[]
    for s in ((M+m)**2*1.001,.004,.007,.02,1.):
        root=np.sqrt(s)
        r=np.sqrt((s-(M+m)**2)*(s-(M-m)**2))/(2*root)
        k=(s-ms*ms)/(2*root)
        x=r*np.sqrt(1-c*c)
        B=np.stack([np.hypot(r,m)+0*c,x,0*c,r*c],axis=-1)
        P=np.stack([np.hypot(r,M)+0*c,-x,0*c,-r*c],axis=-1)
        K=np.array([k,0.,0.,k])
        check=-sum(amplitudes(-P,K,B,M,m,pol)[0] for pol in np.eye(4)[1:3])
        inv=angular(s,c,M,m,ms)
        assert np.min(inv)>0.
        errors.append(float(np.max(abs(check-inv)/inv)))
    assert max(errors)<1e-9,errors
    return max(errors)


def source(g,n,q):
    with gzip.open(GRID/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz','rt') as f:
        a=json.load(f)
    ms=a['checkpoints'][-1]['scalar_pole_over_T']
    M,m=.05,.003
    es=np.hypot(q,ms)
    ep=np.asarray(a['endpoint_parent_momentum_edges'])
    eb=np.asarray(a['endpoint_fermion_momentum_edges'])
    FP=np.asarray(a['endpoint_parent_occupations'])
    FB=np.asarray(a['endpoint_fermion_occupations'])
    h0=(M+m)**2-ms*ms
    cuts=edges_between(np.array([h0/(2*(es+q)),h0/(2*es),h0/(2*(es-q)),.8,1.,1.5,2.,2.5,3.]),.5,4.)
    x,w=roots_legendre(n)
    az,aw=roots_legendre(n+2)
    az,aw=np.pi*(1+az),np.pi*aw
    hk=np.diff(cuts)[:,None]/2
    ks=((cuts[:-1]+cuts[1:])[:,None]/2+hk*x).ravel()
    kw=(hk*w).ravel()
    res=np.zeros(3)
    out=0.
    for k,wk in zip(ks,kw,strict=True):
        ET=k+es
        B=1/np.expm1(k)
        hi=min(1.,(2*k*es-h0)/(2*k*q))
        if hi<=-1:
            continue
        for u,wu in zip((hi-1)/2+(hi+1)*x/2,(hi+1)*w/2,strict=True):
            s=ms*ms+2*k*(es-q*u)
            root=np.sqrt(s)
            r=np.sqrt((s-(M+m)**2)*(s-(M-m)**2))/(2*root)
            P=np.sqrt(k*k+q*q+2*k*q*u)
            gamma=ET/root
            lo=max(np.hypot(.45,m),ET-np.hypot(5.,M))
            high=min(np.hypot(5.,m),ET-np.hypot(.45,M))
            c0=gamma*np.hypot(r,m)
            spread=P/root*r
            u0=max(lo,c0-spread)
            u1=min(high,c0+spread)
            if u1<=u0:
                continue
            para=gamma*(k*(k+q*u)/P-P/ET*k)
            perp=k*q*np.sqrt(1-u*u)/P
            ca,sa=np.array([para,perp])/np.hypot(para,perp)
            ends=np.r_[np.hypot(eb,m),ET-np.hypot(ep,M)]
            for mass in (M,m):
                rho=mass*mass/(r*(np.hypot(r,mass)+r))
                off=rho*np.array([.125,.5,2,8,32,128])
                ends=np.r_[ends,c0+spread*(1-off),c0-spread*(1-off)]
            bins=edges_between(ends,u0,u1)
            h=np.diff(bins)[:,None]/2
            E=((bins[1:]+bins[:-1])[:,None]/2+h*x).ravel()
            ww=(h*w).ravel()
            C=(E-c0)/spread
            cps=C[:,None]*ca+np.sqrt(np.maximum(0.,1-C*C))[:,None]*sa*np.cos(az)
            V=angular(s,np.clip(cps,-1,1),M,m,ms)@aw
            f=occup(np.sqrt((E-m)*(E+m)),eb,FB)
            t=ET-E
            p=occup(np.sqrt((t-M)*(t+M)),ep,FP)
            emission=f*p*(1+B)
            loss=(1-f)*(1-p)*B
            Eqt=1/np.exp(t)/(1+np.exp(-t))
            Eqb=1/np.exp(E)/(1+np.exp(-E))
            # Common T with zero chemical potentials, independently of any
            # nonthermal archive. The reverse/emission ratio is exp(es).
            out=max(out,float(np.max(abs((1-Eqt)*(1-Eqb)*B/(Eqt*Eqb*(1+B))/np.exp(es)-1))))
            pref=k*wk*wu*(r/root)/spread/(256*np.pi**4*es)
            res+=2*pref*np.r_[ww@(V*emission),ww@(V*loss),ww@(V*(emission-loss))]
    assert out<1e-13
    assert abs(res[0]-res[1]-res[2])<1e-10
    return dict(grid=g,q=q,order=n,emission_absorption_and_stimulated_without_couplings=res.tolist(),
                thermal_balance_defect=out,independent_crossed_trace_relative_defect=trace_control(M,m,ms))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--q',type=float,default=4.1e-6)
    p.add_argument('--order',type=int,default=8)
    a=p.parse_args()
    print(json.dumps(source(a.grid,a.order,a.q),indent=2,allow_nan=False))
