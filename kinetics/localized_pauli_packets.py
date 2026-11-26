"""Localized χ-energy logit packet Rayleigh audit for full-Pauli χe contact.

A Gaussian test perturbation of logit(fχ) in x=Eχ/T, projected off
conserved fermionic number, has physical susceptibility S and full 5D
quantum collision Dirichlet D. Rayleigh D/S is a variational upper
bound within the model, not an eigenvalue nor a full continuum gap.

Critical numerical control: narrow IR packets are under-resolved by
ordinary n-point Gauss spanning [0,32T]. Split initial p and k
integrals at 2T to resolve the IR independently from thermal tails.
There is no SM bath dynamics, no gamma_el/CLASS/relic prediction.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import expit
from relativistic_weak_form_2to2 import cm_event,_gauss_span,_input_checks


def packet_rayleigh(center=2.,width=1.,M=.1,T=1.,G=.001,
                    gchi=4.,ge=4.,n=20,cutoff=32.,split=False):
    _input_checks(M,T,T,G,gchi,ge,n,n,n,n,n,cutoff)
    if not all(math.isfinite(z) for z in (center,width)) or width<=0:
        raise ValueError("finite center and strictly positive width")
    p,w=_gauss_span(360,0,70*max(M,T))
    x=np.hypot(p,M)/T
    f=expit(-x)
    weights=gchi/(2*math.pi**2)*w*p*p*f*(1-f)
    def packet(E):
        return np.exp(-.5*((np.asarray(E)/T-center)/width)**2)
    mean=np.dot(weights,packet(np.hypot(p,M)))/np.sum(weights)
    def h(E):
        return packet(E)-mean
    S=float(np.dot(weights,h(np.hypot(p,M))**2))
    if S<=1e-16:
        raise ValueError("packet susceptibility too small; invalid packet")
    if split:
        # 2n-point composite quadrature resolves p~0.2T modes.
        ps=[_gauss_span(n,0,2*T),_gauss_span(n,2*T,cutoff*T)]
        ks=[_gauss_span(n,0,2*T),_gauss_span(n,2*T,cutoff*T)]
        pgrid=np.concatenate([z[0] for z in ps])
        pw=np.concatenate([z[1] for z in ps])
        kgrid=np.concatenate([z[0] for z in ks])
        kw=np.concatenate([z[1] for z in ks])
    else:
        pgrid,pw=_gauss_span(n,0,cutoff*T)
        kgrid,kw=_gauss_span(n,0,cutoff*T)
    mugrid,mw=_gauss_span(n,-1,1)
    ugrid,uw=_gauss_span(n,-1,1)
    angles=2*np.pi*(np.arange(n)+.5)/n
    U=ugrid[:,None];Phi=angles[None,:];wcm=uw[:,None]/n
    D=0.
    for pi,wi in zip(pgrid,pw):
        E=math.hypot(M,pi)
        fchi=float(expit(-E/T))
        h0=float(h(E))
        for k,wk in zip(kgrid,kw):
            fe=float(expit(-k/T))
            for angle,wmu in zip(mugrid,mw):
                s=M*M+2*k*(E-pi*angle)
                flux=1-pi*angle/E
                Eprime,kprime,t,_=cm_event(M,pi,k,angle,U,Phi)
                u=2*M*M-s-t
                delta=s-M*M
                amp=2*G*G*(delta*delta+(u-M*M)**2+2*M*M*t)
                dsdom=amp/(64*math.pi**2*s)
                feq=fchi*fe*(1-expit(-Eprime/T))*(1-expit(-kprime/T))
                dh=h(Eprime)-h0
                D+=wi*pi*pi*wk*k*k*wmu*flux*float(
                    np.sum(wcm*dsdom*feq*dh*dh))
    D*=gchi*ge/(8*math.pi**3)
    return {
      "center_E_over_T":center,"width_E_over_T":width,
      "Rayleigh_GeV":float(D/S),"S_GeV3":S,"D_GeV4":D,
      "n":n,"cutoff":cutoff,"split":split,
      "not_full_spectral_gap":True
    }
