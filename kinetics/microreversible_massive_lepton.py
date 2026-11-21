"""Directly reciprocal full-Pauli chi + MASSIVE charged lepton vector 2->2.

Spin-averaged |M|²=2G²[(s-M²-m²)²+(u-M²-m²)²+2t(M²+m²)].
Exact invariant q-transfer kernel gchi ge/(256 pi^5) ∫dq <|M|²>.
No posthoc rate symmetrization. Electron+positron mass and pair-threshold
phase space included, but not QED plasma or e+e- annihilation/charge
neutrality, UV/EW gauge completion, FRW, Bose 1<->2 or 86-eV relic.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit
from microreversible_chie_lattice import sources

ELECTRON_MASS_GEV=0.00051099895


def massive_vector_kernel(M,m,G,E,eps,Ep,epsp,gchi=4.,ge=4.,nq=12):
    v=(M,m,G,E,eps,gchi,ge)
    if (not all(math.isfinite(float(z)) for z in v) or M<=0 or m<0
        or G<0 or E<=M or eps<=m or gchi<=0 or ge<=0
        or type(nq) is not int or nq<2):
        raise ValueError("physical initial states, masses, coupling required")
    Ep=np.asarray(Ep,float);epsp=np.asarray(epsp,float)
    if (Ep.shape!=epsp.shape or np.any(~np.isfinite(Ep)) or
        np.any(~np.isfinite(epsp)) or
        np.any(np.abs(E+eps-Ep-epsp)>1e-10*max(1.,E+eps))):
        raise ValueError("outgoing physical energy conservation required")
    p=math.sqrt(E*E-M*M);k=math.sqrt(eps*eps-m*m)
    pp=np.sqrt(np.maximum(0,Ep*Ep-M*M))
    kp=np.sqrt(np.maximum(0,epsp*epsp-m*m))
    low=np.maximum(np.abs(p-pp),np.abs(k-kp))
    high=np.minimum(p+pp,k+kp)
    ok=(high>low)&(Ep>M)&(epsp>m)
    low=np.where(ok,low,0.);high=np.where(ok,high,0.)
    x,w=leggauss(nq)
    q=low[...,None]+(high-low)[...,None]*(1+x)/2
    qq=np.where(q>0,q,1.)
    a=np.clip((p*p+qq*qq-pp[...,None]**2)/(2*p*qq),-1,1)
    b=np.clip((kp[...,None]**2-k*k-qq*qq)/(2*k*qq),-1,1)
    D0=2*(E*eps-p*k*a*b)
    t=(E-Ep[...,None])**2-qq*qq
    D2=D0*D0+2*p*p*k*k*(1-a*a)*(1-b*b)
    amp=2*G*G*(2*D2+2*D0*t+t*t+2*(M*M+m*m)*t)
    if np.any(amp[ok]<-1e-12):
        raise ArithmeticError("unphysical negative squared amplitude")
    return gchi*ge/(256*math.pi**5)*np.sum(
        np.maximum(0.,amp)*w,axis=-1)*(high-low)/2*ok


def massive_chie_lattice(*,M=.01,m=ELECTRON_MASS_GEV,T=.0001,G=.001,
                         dE=.000025,nchi=100,ne=100,
                         gchi=4.,ge=4.,nq=10):
    v=(M,m,T,G,dE,gchi,ge)
    if (not all(math.isfinite(float(z)) for z in v) or
        M<=0 or m<0 or T<=0 or G<0 or dE<=0 or gchi<=0 or ge<=0
        or any(type(n) is not int or n<3 for n in (nchi,ne))
        or nchi*ne>45000):
        raise ValueError("physical finite-energy lattice required")
    E=M+(np.arange(nchi)+.5)*dE
    K=m+(np.arange(ne)+.5)*dE
    pchi=np.sqrt(np.maximum(0,(M+dE*np.arange(nchi+1))**2-M*M))
    pe=np.sqrt(np.maximum(0,(m+dE*np.arange(ne+1))**2-m*m))
    Vchi=gchi/(6*math.pi**2)*np.diff(pchi**3)
    Ve=ge/(6*math.pi**2)*np.diff(pe**3)
    ii=[];aa=[];jj=[];bb=[];ww=[]
    for i,e in enumerate(E):
        for a,k in enumerate(K):
            j=np.arange(i+1,min(nchi,i+a+1),dtype=int)
            if not len(j):continue
            b=a+i-j
            rates=massive_vector_kernel(M,m,G,float(e),float(k),
                 E[j],K[b],gchi=gchi,ge=ge,nq=nq)*dE**3
            good=rates>0
            ii.extend([i]*int(np.sum(good)))
            aa.extend([a]*int(np.sum(good)))
            jj.extend(j[good].tolist())
            bb.extend(b[good].tolist())
            ww.extend(rates[good].tolist())
    i=np.array(ii,int);a=np.array(aa,int)
    j=np.array(jj,int);b=np.array(bb,int)
    W=np.array(ww,float)
    if not np.allclose(E[i]+K[a],E[j]+K[b],atol=1e-12,rtol=1e-13):
        raise ArithmeticError("quartet energy mismatch")
    return dict(Echi=E,Ee=K,Vchi=Vchi,Ve=Ve,i=i,a=a,j=j,b=b,W=W,
        params=dict(M=M,me=m,T=T,G=G,dE=dE,nchi=nchi,ne=ne,
                    gchi=gchi,ge=ge,nq=nq),
        direct_microreversible_unequal_mass_amplitude=True,
        no_baryonic_electrons_or_plasma_annihilation=True)


def thermal_heat_example(*,T=.0001,m=ELECTRON_MASS_GEV,
                         M=.01,G=.001,n=100,nq=8):
    g=massive_chie_lattice(M=M,m=m,T=T,G=G,dE=T/4,
                           nchi=n,ne=n,nq=nq)
    f=expit(((M-2*T)-g["Echi"])/(.8*T))
    e=expit(-g["Ee"]/T)
    r=sources(g,f,e)
    return dict(T=T,me=m,Mchi=M,G=G,Qchi=r["Qchi_GeV5"],
        Qe=r["Qe_GeV5"],entropy=r["entropy_GeV4"],
        energy_residual=r["energy_residual_GeV5"],
        scope="ideal charge-symmetric e+/e- FD bath, mu_e=0, mu_chi=M-2T; no annihilation")
