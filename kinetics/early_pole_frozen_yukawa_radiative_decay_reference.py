"""Frozen real radiative 1↔3 contribution on strictly positive photon cuts."""
import argparse
import gzip
import json

import numpy as np
from scipy.special import roots_legendre

from early_pole_frozen_yukawa_photoproduction_source_reference import edges_between,parent_input
from early_pole_hard_yukawa_photoproduction_reference import amplitudes,dot,invariant_form
from early_pole_native_aperture_witness_reference import GRID,occup


def trace_control(M,m,ms):
    P=np.array([M,0.,0.,0.])
    defects=[]
    for frac in (.02,.2,.5,.8,.98):
        s=(m+ms)**2+frac*(M*M-(m+ms)**2)
        root=np.sqrt(s)
        r=np.sqrt((s-(m+ms)**2)*(s-(m-ms)**2))/(2*root)
        k=(M*M-s)/(2*M)
        E=np.hypot(r,m)
        beta=-k/(M-k)
        gamma=(M-k)/root
        c=np.linspace(-.97,.97,9)
        B=np.stack([gamma*(E+beta*r*c),r*np.sqrt(1-c*c),0*c,gamma*(r*c+beta*E)],axis=-1)
        K=np.array([k,0.,0.,k])
        x=invariant_form(-M*k,-dot(B,K),dot(P,B),M,m)
        y=sum(amplitudes(P,-K,B,M,m,v)[0] for v in np.eye(4)[1:3])
        assert np.min(x)>0
        z=P-K-B
        assert np.max(abs(dot(z,z)-ms*ms))<1e-17
        defects.append(float(np.max(abs(x-y)/x)))
    assert max(defects)<1e-9
    return max(defects)


def source(g,n,q,witness_mix=0.):
    with gzip.open(GRID/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz','rt') as f:
        a=json.load(f)
    ms=a['checkpoints'][-1]['scalar_pole_over_T']
    M,m=.05,.003
    es=np.hypot(ms,q)
    pe,FP=parent_input(g,a,witness_mix)
    eb=np.asarray(a['endpoint_fermion_momentum_edges'])
    FB=np.asarray(a['endpoint_fermion_occupations'])
    edges=edges_between(eb,.45,5.)
    ue=np.hypot(pe,M)
    x,w=roots_legendre(n)
    az,aw=roots_legendre(n+2)
    az,aw=np.pi*(1+az),np.pi*aw
    hk=np.diff(edges)[:,None]/2
    pp=((edges[:-1]+edges[1:])[:,None]/2+hk*x).ravel()
    wp=(hk*w).ravel()
    res=np.zeros(3)
    terr=0.
    for p,W in zip(pp,wp,strict=True):
        E= np.hypot(p,m)
        f=occup(np.array(p),eb,FB)
        A=E+es
        lower=max(-1.,(ms*ms+m*m+2*E*es-M*M)/(2*p*q))
        if lower>=1.:
            continue
        for u,Wu in zip((1+lower)/2+(1-lower)*x/2,(1-lower)*w/2,strict=True):
            s=m*m+ms*ms+2*(E*es-p*q*u)
            Z=np.sqrt(p*p+q*q+2*p*q*u)
            half=(M*M-s)/2
            if half<=0.:
                continue
            lo=max(.5,half/(A+Z),np.hypot(.45,M)-A)
            hi=min(4.,half*(A+Z)/s,np.hypot(5.,M)-A)
            if hi<=lo:
                continue
            rot,rad=np.array([p+q*u,q*np.sqrt(1-u*u)])/Z
            bins=edges_between(ue-A,lo,hi)
            h=np.diff(bins)[:,None]/2
            k=((bins[1:]+bins[:-1])[:,None]/2+h*x).ravel()
            wk=(h*w).ravel()
            c=(A-half/k)/Z
            cps=c[:,None]*rot+np.sqrt(np.maximum(0.,1-c*c))[:,None]*rad*np.cos(az)
            bk=k[:,None]*(m*m/(E+p)+p*np.maximum(0.,1-cps))
            V=invariant_form(-half,-bk,m*m+(s-m*m-ms*ms)/2+bk,M,m)@aw
            Et=A+k
            F=occup(np.sqrt((Et-M)*(Et+M)),pe,FP)
            B=1/np.expm1(k)
            emit=F*(1-f)*(1+B)
            loss=f*(1-F)*B
            ppq=1/(1+np.exp(Et))
            fbq=1/(1+np.exp(E))
            terr=max(terr,float(np.max(abs((fbq*(1-ppq)*B)/(ppq*(1-fbq)*(1+B))/np.exp(es)-1))))
            measure=p*p*W*Wu/(256*np.pi**4*E*es*Z)*wk
            res+=2*np.r_[measure@(V*emit),measure@(V*loss),measure@(V*(emit-loss))]
    assert terr<2e-13
    assert abs(res[0]-res[1]-res[2])<1e-10
    return dict(grid=g,q=q,order=n,witness_mix=witness_mix,
                emission_absorption_and_stimulated_without_couplings=res.tolist(),
                thermal_balance_defect=terr,independent_trace_relative_defect=trace_control(M,m,ms))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--q',type=float,default=4.1e-6)
    p.add_argument('--order',type=int,default=8)
    p.add_argument('--witness-mix',type=float,default=0.)
    a=p.parse_args()
    print(json.dumps(source(a.grid,a.order,a.q,a.witness_mix),indent=2,allow_nan=False))
