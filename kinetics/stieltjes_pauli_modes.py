"""Orthogonal Stieltjes recurrence for up to 11 nonconstant Pauli Ritz modes.

FD susceptibility measure x=E/T, weight gchi*p^2*f*(1-f)/(2pi^2).
Polynomials are evaluated by stable three-term recurrence instead of
high-degree monomial subtraction. Full 5D Pauli chi-e Dirichlet form
evaluated independently from previous 2-5-mode polynomial implementation.

Finite Galerkin Ritz eigenvalues can be upper bounds on variational slow
eigenmodes; they DO NOT establish a full Boltzmann spectral gap, a real
kinetic gamma_el, finite SM bath or physical 86-eV relic.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import expit
from scipy.linalg import eigh
from relativistic_weak_form_2to2 import cm_event,_gauss_span,_input_checks


def stieltjes_basis(M=.1,T=1.,mu=0.,gchi=4.,modes=8):
    if (type(modes) is not int or not 1<=modes<=12 or
        M<=0 or T<=0 or gchi<=0 or
        not all(math.isfinite(z) for z in (M,T,mu,gchi))):
        raise ValueError("positive finite thermal parameters and 1..12 modes")
    p,pw=_gauss_span(450,0,90*max(M,T))
    x=np.hypot(p,M)/T
    f=expit((mu-T*x)/T)
    weight=pw*p**2*f*(1-f)*gchi/(2*math.pi**2)
    norm0=math.sqrt(np.dot(weight,np.ones_like(x)))
    qprev=np.zeros_like(x)
    q=np.ones_like(x)/norm0
    alpha=[];beta=[];prev_beta=0.
    for _ in range(modes):
        a=np.dot(weight*x*q,q)
        residual=(x-a)*q-prev_beta*qprev
        b=math.sqrt(np.dot(weight,residual**2))
        alpha.append(a);beta.append(b)
        if b<=1e-14:raise ArithmeticError("Degenerate Stieltjes polynomial")
        qprev,q=q,residual/b
        prev_beta=b
    def evaluate(E):
        x=np.asarray(E)/T
        prev=np.zeros_like(x)
        curr=np.ones_like(x)/norm0
        out=[]
        for n in range(modes):
            out.append(curr)
            next_=((x-alpha[n])*curr-(beta[n-1]*prev if n else 0))/beta[n]
            prev,curr=curr,next_
        return np.array(out[1:])
    return evaluate,dict(alpha=alpha,beta=beta,constant_norm=norm0),p,weight


def pauli_ritz(M=.1,T=1.,mu=0.,G=.001,gchi=4.,ge=4.,
                nonconstant_modes=7,n=20,cutoff=32.):
    if type(nonconstant_modes) is not int or not 1<=nonconstant_modes<=11:
        raise ValueError("nonconstant mode count 1..11")
    _input_checks(M,T,T,G,gchi,ge,n,n,n,n,n,cutoff)
    evalm,params,p_therm,weight=stieltjes_basis(
         M,T,mu,gchi,nonconstant_modes+1)
    pgrid,pw=_gauss_span(n,0,cutoff*T)
    kgrid,kw=_gauss_span(n,0,cutoff*T)
    mu_grid,mw=_gauss_span(n,-1,1)
    u_grid,uw=_gauss_span(n,-1,1)
    phi=2*np.pi*(np.arange(n)+.5)/n
    uu=u_grid[:,None];ph=phi[None,:];wcm=uw[:,None]/n
    D=np.zeros((nonconstant_modes,nonconstant_modes))
    for p,wp in zip(pgrid,pw):
        E=math.hypot(M,p)
        fi=float(expit((mu-E)/T))
        hin=evalm(E)
        for k,wk in zip(kgrid,kw):
            fe=float(expit(-k/T))
            for cost,wmu in zip(mu_grid,mw):
                s=M*M+2*k*(E-p*cost)
                flux=1-p*cost/E
                Ep,kp,t,_=cm_event(M,p,k,cost,uu,ph)
                umand=2*M*M-s-t
                delta=s-M*M
                amp=2*G*G*(delta**2+(umand-M*M)**2+2*M*M*t)
                dsdom=amp/(64*math.pi**2*s)
                fp=expit((mu-Ep)/T)
                gp=expit(-kp/T)
                Feq=fi*fe*(1-fp)*(1-gp)
                dh=evalm(Ep)-hin[:,None,None]
                weights=wp*p*p*wk*k*k*wmu*flux*wcm*dsdom*Feq
                for a in range(nonconstant_modes):
                    for b in range(a,nonconstant_modes):
                        D[a,b]+=float(np.sum(weights*dh[a]*dh[b]))
    D=D+np.triu(D,1).T
    D*=gchi*ge/(8*math.pi**3)
    therm_basis=evalm(np.hypot(p_therm,M))
    S=np.einsum("p,ip,jp->ij",weight,therm_basis,therm_basis)
    eig=eigh(D,S,eigvals_only=True)
    return {
      "nonconstant_modes":nonconstant_modes,"rates_GeV":eig.tolist(),
      "basis_orthogonality_max":float(
          np.max(np.abs(S-np.eye(nonconstant_modes)))),
      "cutoff":cutoff,"n":n,
      "rate_means_no_spectral_gap_proof":True
    }


if __name__=="__main__":
    import json
    print(json.dumps(pauli_ritz(nonconstant_modes=11,n=24),indent=2))
