"""Two-mode Galerkin linearized full-Pauli chi+e quantum 2-to-2 collision.

At equilibrium fchi=FD((E-mu)/T), fe=FD(k/T), a logit shift
h(E)=delta log[f/(1-f)] generates the linear response. The symmetric
Dirichlet form is
  D_ab=(gchi*ge/2) int[d3p d3k/(2pi)^6] vM
          int dOmega_cm dsigma/dOmega_cm F_eq [h_a(E')-h_a(E)]
                                                   [h_b(E')-h_b(E)].
The equilibrium susceptibility is
  S_ab=gchi/(2pi^2) int dp p^2 f(1-f) h_a(E)h_b(E).
The 2 generalized positive eigenvalues Dv=lambda Sv are only TWO
Ritz values of the restricted number-conserving isotropic logit
perturbation space. Not the full linearized Boltzmann spectrum, not
a real phenomenological BGK gamma, no finite SM feedback, no relic.

Electron massless, e± fixed ideal FD, spin-averaged vector contact
in broken-EW EFT, G[GeV^-2], M,T[GeV], lambda[GeV].
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import expit
from scipy.integrate import quad
from scipy.linalg import eigh
from relativistic_weak_form_2to2 import cm_event,_gauss_span,_input_checks


def fd_logit_basis(M,T,mu=0.):
    if M<=0 or T<=0 or not all(math.isfinite(v) for v in (M,T,mu)):
        raise ValueError("positive M,T and finite chemical potential required")
    def moment(j):
        def calc(p):
            E=math.hypot(p,M)
            f=float(expit((mu-E)/T))
            return p*p*f*(1-f)*(E/T)**j
        return quad(calc,0,60*max(M,T),epsabs=1e-11,epsrel=1e-10,
                    limit=160)[0]
    mm=[moment(i) for i in range(5)]
    mean=mm[1]/mm[0]
    var=mm[2]/mm[0]-mean*mean
    if var<=0: raise ArithmeticError("Nonpositive finite FD thermal variance")
    central3=mm[3]/mm[0]-3*mean*mm[2]/mm[0]+2*mean**3
    skewcoef=central3/var
    def h1(E):
        x=np.asarray(E)/T-mean
        return x
    def h2(E):
        x=np.asarray(E)/T-mean
        return x*x-var-skewcoef*x
    norm1=var*mm[0]
    norm2=mm[4]-4*mean*mm[3]+6*mean*mean*mm[2]-4*mean**3*mm[1]+mean**4*mm[0]
    norm2-=mm[0]*var**2+mm[0]*central3**2/var
    if norm2<=0: raise ArithmeticError("Nonpositive second mode norm")
    return h1,h2,(norm1,norm2),{
      "mean_E_over_T":mean,"variance":var,"skewcoefficient":skewcoef}


def linearized_modes(M=.1,T=1.,mu=0.,G=.001,gchi=4.,ge=4.,
                     n=12,cutoff=24.):
    _input_checks(M,T,T,G,gchi,ge,n,n,n,n,n,cutoff)
    h1,h2,norm,details=fd_logit_basis(M,T,mu)
    # Independent continuum FD susceptibility (integrates up to 60 max[M,T]).
    S=(gchi/(2*math.pi**2))*np.diag(norm)
    D=np.zeros((2,2),float)
    pgrid,pw=_gauss_span(n,0,cutoff*T)
    kgrid,kw=_gauss_span(n,0,cutoff*T)
    mugrid,mw=_gauss_span(n,-1,1)
    ugrid,uw=_gauss_span(n,-1,1)
    phigrid=2*np.pi*(np.arange(n)+.5)/n
    uu=ugrid[:,None]; ph=phigrid[None,:]; wcm=uw[:,None]/n
    pref=gchi*ge/(8*math.pi**3)
    for p,wp in zip(pgrid,pw):
        E=math.hypot(M,p)
        fi=float(expit((mu-E)/T))
        h=np.array([float(h1(E)),float(h2(E))])
        for k,wk in zip(kgrid,kw):
            gi=float(expit(-k/T))
            for cosi,wmu in zip(mugrid,mw):
                s=M*M+2*k*(E-p*cosi)
                flux=1-p*cosi/E
                Eprime,eeprime,t,_=cm_event(M,p,k,cosi,uu,ph)
                A=s-M*M
                us=2*M*M-s-t
                amp=2*G*G*(A*A+(us-M*M)**2+2*M*M*t)
                dsdom=amp/(64*math.pi**2*s)
                fprime=expit((mu-Eprime)/T)
                eprime=expit(-eeprime/T)
                Feq=fi*gi*(1-fprime)*(1-eprime)
                dh=np.array([h1(Eprime)-h[0],h2(Eprime)-h[1]])
                common=(wp*p*p*wk*k*k*wmu*flux)*wcm*dsdom*Feq
                for i in range(2):
                    for j in range(i,2):
                        D[i,j]+=float(np.sum(common*dh[i]*dh[j]))
    D[1,0]=D[0,1]
    D*=pref
    eig,eigvec=eigh(D,S)
    return {
      "D_matrix_GeV4":D.tolist(),"M_matrix_GeV3":S.tolist(),
      "generalized_rates_GeV":eig.tolist(),
      "normalized_eigenvectors":eigvec.tolist(),
      "thermal_energy_mode_rate_GeV":float(D[0,0]/S[0,0]),
      "coupling_offdiagonal_fraction":
           float(D[0,1]/math.sqrt(D[0,0]*D[1,1])) if D[0,0]*D[1,1]>0 else 0.,
      "basis":details,
      "parameter_scope":dict(M=M,T=T,mu=mu,G=G,gchi=gchi,ge=ge,
                             n=n,cutoff=cutoff),
      "not_full_spectrum_nor_kinetic_evolution":True
    }


if __name__=="__main__":
    import json
    print(json.dumps(linearized_modes(n=16,cutoff=28.),indent=2))
