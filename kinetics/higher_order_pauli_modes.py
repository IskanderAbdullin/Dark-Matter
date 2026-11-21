"""FD-weighted orthogonal-polynomial Galerkin Ritz hierarchy (2..5 modes).

All basis functions are orthogonal to the conserved χ number direction.
More modes can only lower the smallest Ritz eigenvalue in exact arithmetic.
Two/three/four/five mode eigenvalues are NOT the full nonthermal kinetic
spectrum, and do not establish a nonzero true spectral gap or a BGK gamma.
Electron-only vector-contact full-Pauli 5D collision, fixed electron FD.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import expit
from scipy.linalg import eigh
from relativistic_weak_form_2to2 import cm_event,_gauss_span,_input_checks


def number_orthogonal_polynomial_basis(M,T,mu=0.,gchi=4.,modes=3):
    if (type(modes) is not int or not 1<=modes<=5 or
        not (M>0 and T>0 and gchi>0) or
        not all(math.isfinite(x) for x in (M,T,mu,gchi))):
        raise ValueError("finite M,T,g>0 and number of modes 1..5 required")
    pp,pw=_gauss_span(280,0,60*max(M,T))
    E=np.hypot(pp,M);x=E/T
    f=expit((mu-E)/T)
    weights=pw*pp**2*f*(1-f)*gchi/(2*math.pi**2)
    coeffs=[]
    for degree in range(modes+1):
        c=np.zeros(modes+1);c[degree]=1.
        # Twice-reorthogonalized Gram-Schmidt for nonnegative FD measure.
        for _ in range(2):
            for v in coeffs:
                prod=float(np.dot(weights,
                   np.polynomial.polynomial.polyval(x,c)*
                   np.polynomial.polynomial.polyval(x,v)))
                c-=prod*v
        length=math.sqrt(float(np.dot(weights,
                         np.polynomial.polynomial.polyval(x,c)**2)))
        if length<=1e-10:
            raise ArithmeticError("degenerate high-order thermal mode")
        coeffs.append(c/length)
    # Exclude constant null direction. Check charge and mutual overlaps.
    return np.array(coeffs[1:]),np.array(coeffs[0]),weights,pp


def ritz_hierarchy(M=.1,T=1.,mu=0.,G=.001,gchi=4.,ge=4.,
                   modes=3,n=16,cutoff=28.):
    _input_checks(M,T,T,G,gchi,ge,n,n,n,n,n,cutoff)
    basis,const,weight,pp=number_orthogonal_polynomial_basis(
       M,T,mu,gchi,modes)
    dim=modes
    D=np.zeros((dim,dim))
    pgrid,pw=_gauss_span(n,0,cutoff*T)
    kgrid,kw=_gauss_span(n,0,cutoff*T)
    mugrid,mw=_gauss_span(n,-1,1)
    ugrid,uw=_gauss_span(n,-1,1)
    phi=2*np.pi*(np.arange(n)+.5)/n
    u=ugrid[:,None];ph=phi[None,:];wcm=uw[:,None]/n
    pref=gchi*ge/(8*math.pi**3)
    def eval_modes(E):
        return np.asarray([
          np.polynomial.polynomial.polyval(np.asarray(E)/T,c)
          for c in basis])
    for p,wp in zip(pgrid,pw):
        E=math.hypot(M,p)
        fchi=float(expit((mu-E)/T))
        hin=eval_modes(E)
        for k,wk in zip(kgrid,kw):
            fe=float(expit(-k/T))
            for cost,wmu in zip(mugrid,mw):
                s=M*M+2*k*(E-p*cost)
                flux=1-p*cost/E
                Ep,kp,t,_=cm_event(M,p,k,cost,u,ph)
                delta=s-M*M;umat=2*M*M-s-t
                amp=2*G*G*(delta*delta+(umat-M*M)**2+2*M*M*t)
                dsdom=amp/(64*math.pi**2*s)
                fp=expit((mu-Ep)/T);gp=expit(-kp/T)
                Feq=fchi*fe*(1-fp)*(1-gp)
                dh=eval_modes(Ep)-hin[:,None,None]
                common=(wp*p*p*wk*k*k*wmu*flux)*wcm*dsdom*Feq
                for a in range(dim):
                    for b in range(a,dim):
                        D[a,b]+=float(np.sum(common*dh[a]*dh[b]))
    D=D+np.triu(D,1).T
    D*=pref
    # Susceptibility independently evaluated on dense finite FD support.
    qx=np.hypot(pp,M)/T
    S=np.array([
      [float(np.dot(weight,
       np.polynomial.polynomial.polyval(qx,a)*
       np.polynomial.polynomial.polyval(qx,b))) for b in basis]
      for a in basis])
    eig=eigh(D,S,eigvals_only=True)
    null=float(max(abs(np.dot(weight,
       np.polynomial.polynomial.polyval(qx,const)*
       np.polynomial.polynomial.polyval(qx,b))) for b in basis))
    return {
       "modes":modes,"Ritz_rates_GeV":eig.tolist(),
       "D_GeV4":D.tolist(),"S_GeV3":S.tolist(),
       "number_orthogonality_max":null,
       "basis_orthonormality_max":float(np.max(np.abs(S-np.eye(dim)))),
       "n":n,"cutoff":cutoff,
       "warning":"finite Galerkin Ritz rates, not full spectral gap, no time evolution or DM relic"
    }


if __name__=="__main__":
    import json
    print(json.dumps(ritz_hierarchy(n=20,cutoff=32,modes=5),indent=2))
