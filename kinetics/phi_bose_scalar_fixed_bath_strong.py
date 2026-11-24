"""Conditional on-shell Yukawa scalar mode strong equation in fixed fermion bath.

Physical scalar phase-space DOS and inverse parent kinematic support are used,
yielding exact positive per-mode semigroup. The two fermion spectra are HELD
FIXED: this is not a self-consistent three-species strong operator and makes
no SM/FRW, relic, CLASS or condensate prediction.
"""
import math
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.special import expit
from kinetics.phi_bose_physical_cutcell_weak import (
    M1, MPHI, G1, KSTAR, ESTAR, WIDTH, shell, scalar_edges,
)


def threshold(k,pmax):
    """Smallest physical parent momentum supporting a scalar with energy k."""
    lower,upper=shell(pmax)[1:3]
    if k<lower or k>upper:return None
    if k==ESTAR:return 0.
    edge=1 if k<ESTAR else 2
    return brentq(lambda p:shell(p)[edge]-k,0.,pmax,xtol=2e-13)


def mode_coefficients(k,pmax,nparent=24,T1=.9,T2=1.,mu1=0.,mu2=0.,gphi=1.):
    """(S,D,S1,D1): db/dt=S-D*b at fixed FD1/FD2.

    S and D in GeV, E1-weighted S1 and D1 in GeV². Includes division by
    *physical* real-scalar energy DOS gphi*k*sqrt(k²-mphi²)/(2π²).
    Parent measure uses physical g1 Γrest M1/E1, exact on-shell angle cut.
    """
    assert k>MPHI and nparent>=2 and pmax>0 and T1>0 and T2>0 and gphi>0
    p0=threshold(k,pmax)
    if p0 is None:return 0.,0.,0.,0.
    x,w=leggauss(nparent)
    p=p0+(pmax-p0)*(x+1)/2
    dp=(pmax-p0)*w/2
    e1=np.hypot(p,M1)
    e2=e1-k
    f1=expit(-(e1-mu1)/T1)
    f2=expit(-(e2-mu2)/T2)
    weight=dp*p/e1
    production=f1*(1-f2)
    depletion=f2-f1
    pphi=math.sqrt((k-MPHI)*(k+MPHI))
    pref=G1*WIDTH*M1*M1/(2*gphi*KSTAR*k*pphi)
    return tuple(float(pref*np.dot(weight,z)) for z in (
        production,depletion,e1*production,e1*depletion))


def exact_mode_step(k,pmax,dt,b0,nparent=24,T1=.9,T2=1.,mu1=0.,mu2=0.):
    """Exact nonnegative scalar occupation for fixed fermions and finite dt>=0.

    For D>=0, b=b0 exp(-Ddt)+S(1-exp(-Ddt))/D. For D<0 the
    stimulated solution still stays positive (but may grow/overflow).
    """
    assert dt>=0 and b0>=0
    S,D,_,_=mode_coefficients(k,pmax,nparent,T1,T2,mu1,mu2)
    x=-D*dt
    if x>700:raise OverflowError("Positive Bose stimulation exceeds float range")
    exprel=math.expm1(x)/x if x!=0 else 1.
    return b0*math.exp(x)+S*dt*exprel


def integrate_weak(pmax,ncells=16,nk=32,nparent=32,
                   T1=.9,T2=1.,mu1=0.,mu2=0.,muphi=0.):
    """Independent integration-order-reversed Rphi,Qphi,Q1,Q2 moments."""
    assert muphi<MPHI
    edges=scalar_edges(pmax,ncells)
    kmin,kmax=shell(pmax)[1:3]
    x,w=leggauss(nk)
    bycell=np.zeros((ncells,4))
    for j,(left,right) in enumerate(zip(edges[:-1],edges[1:])):
        a=max(left,kmin);b=min(right,kmax)
        if b<=a:continue
        cuts=[a,ESTAR,b] if a<ESTAR<b else [a,b]
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            ln=math.log(hi/lo)
            energies=lo*np.exp((x+1)*ln/2)
            weights=w*energies*ln/2
            for k,dk in zip(energies,weights):
                S,D,S1,D1=mode_coefficients(float(k),pmax,nparent,T1,T2,mu1,mu2)
                bose=1/math.expm1((k-muphi)/T2)
                r=S-D*bose
                p1=S1-D1*bose
                dos=k*math.sqrt((k-MPHI)*(k+MPHI))/(2*math.pi**2)
                weight=dk*dos
                bycell[j]+=(weight*r,weight*k*r,
                             -weight*p1,weight*(p1-k*r))
    total=np.sum(bycell,axis=0)
    return dict(Rphi=float(total[0]),Qphi=float(total[1]),
                Q1=float(total[2]),Q2=float(total[3]),
                energy_defect=float(sum(total[1:])),bycell=bycell)


def firstcell_exact_delta_moments(n,pmax,dt,nk=48,nparent=32):
    """Actual physical Bose first-cell two-moment delta, fixed FD baths.

    Independently integrates b(k,t)-b(k,0), without linearization in dt.
    Outside finite parent support, mode occupation is exactly unchanged.
    """
    upper=scalar_edges(pmax,n)[1]
    low=shell(pmax)[1]
    high=min(upper,shell(pmax)[2])
    if high<=low:return np.zeros(2)
    x,w=leggauss(nk)
    ln=math.log(high/low)
    ks=low*np.exp((x+1)*ln/2)
    dks=w*ks*ln/2
    M=np.zeros(2)
    for k,dk in zip(ks,dks):
        S,D,_,_=mode_coefficients(float(k),pmax,nparent)
        b0=1/math.expm1(k)
        x=-D*dt
        exprel=math.expm1(x)/x if x!=0 else 1.
        db=dt*(S-D*b0)*exprel
        dos=k*math.sqrt((k-MPHI)*(k+MPHI))/(2*math.pi**2)
        hR=(k-MPHI)/(upper-MPHI)
        M+=dk*dos*db*np.array([1-hR,hR])
    return M


def audit():
    from kinetics.phi_bose_physical_cutcell_weak import analytic_weak
    rows=[]
    for P in (2.,8.,16.):
        ref=analytic_weak(P)
        for nr in (8,16,32,48):
            got=integrate_weak(P,16,nr,nr)
            err={key:abs(got[key]/ref[key]-1)
                 for key in ("Rphi","Qphi","Q1","Q2")}
            rows.append(dict(P=P,nk=nr,nparent=nr,errors=err,
                             Rphi=got["Rphi"],Qphi=got["Qphi"],
                             Q1=got["Q1"],Q2=got["Q2"],
                             firstcell_Rphi=float(got["bycell"][0,0])))
    return dict(schema_version=1,
                scope="conditional fixed FD1/FD2 scalar mode semigroup, not self-consistent strong three-species kinetics",
                rows=rows,class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
