"""Exact tree-level vacuum QED e-e+ -> gamma gamma controls and FD LOSS rate.

The thermal rate integrates initial-state FD e± only: NO final-state
stimulated photon Bose factors, NO inverse gamma gamma -> e-e+, NO
chemical-relaxation or full quantum Boltzmann dynamics. This is a
deliberately incomplete collision benchmark, NOT an annihilation freezeout.
Natural units, energies GeV, sigma GeV^-2, reaction density GeV^4.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit

ALPHA=1/137.035999084


def sigma_v_cm(s, m=1., alpha=ALPHA):
    """Exact spin-averaged vacuum 2γ sigma*v_Moller in the pair COM.

    v_Moller=2β; includes identical-photon 1/2. For β<1e-3
    use independently expanded threshold series to prevent cancellation.
    """
    s=np.asarray(s,dtype=float)
    if m<=0 or alpha<=0 or np.any(~np.isfinite(s)) or np.any(s<4*m*m):
        raise ValueError("finite s>=4m², positive m and alpha required")
    b2=np.maximum(0.,1.-4*m*m/s)
    b=np.sqrt(b2)
    bs=np.where(b<.001,1.,b)
    L=2*np.arctanh(np.minimum(bs,1.-1e-15))
    bracket=((3.-b2*b2)/bs)*L-2*(2-b2)
    exact=(1-b2)*bracket/2
    threshold=1+b2-2.4*b2*b2+(52/105)*b2**3+(4/105)*b2**4
    return np.pi*alpha**2/m**2*np.where(b<.001,threshold,exact)


def sigma_v_cm_angular(s,m=1.,alpha=ALPHA,n=128):
    """Independent CM spin-pol amplitude angular integration (no σ formula).

    |M|²_av=4e⁴[(2−D)/D+2(1−β²)/D−2(1−β²)²/D²],
    D=1−β²cos²θ. Factor 1/2 for identical final photons included.
    """
    beta=math.sqrt(1-4*m*m/s)
    c,w=leggauss(n)
    D=1-beta*beta*c*c
    H=(2-D)/D+2*(1-beta*beta)/D-2*(1-beta*beta)**2/D**2
    return 2*math.pi*alpha**2/s*float(np.dot(w,H))


def sigma_v_geant_rest_reference(s,m=1.,alpha=ALPHA):
    """Independent Heitler / Geant4 electron-at-rest cross section.

    Converted back to Lorentz-invariant σ(s)*CM Moller v. Use moderate
    β>0.01 to avoid subtractive threshold cancellation.
    """
    gamma=s/(2*m*m)-1
    if gamma<=1: raise ValueError("nonthreshold reference only")
    rt=math.sqrt(gamma*gamma-1)
    sigma=math.pi*(alpha/m)**2/(gamma+1)*(
       (gamma*gamma+4*gamma+1)/(gamma*gamma-1)*math.acosh(gamma)
       -(gamma+3)/rt)
    beta=math.sqrt(1-4*m*m/s)
    return sigma*2*beta


def initial_fd_loss_density(T,m=1.,alpha=ALPHA,*,n=32,cut=32):
    """Initial FD e−e+->2γ loss R per volume/time (units energy^4).

    Exact vacuum QED σv_M combined with FD occupations, WITHOUT
    final-state Bose stimulation or inverse photon pair production.
    This object is NOT a net plasma collision operator.
    """
    if not(0<T and m>0 and n>=8 and cut>0): raise ValueError
    x,wx=leggauss(n)
    q=cut*(x+1)/2
    W=cut*wx/2
    mu,wmu=leggauss(n)
    E=np.sqrt(m*m+(T*q)**2)
    f=expit(-E/T)
    accum=0.
    for i,qi in enumerate(q):
        for j,qj in enumerate(q):
            s=2*m*m+2*(E[i]*E[j]-T*T*qi*qj*mu)
            # σ v_M(lab)= [σ v_M(COM)]*s/(4 E_i E_j).
            sigma_v=sigma_v_cm(s,m,alpha)*s/(4*E[i]*E[j])
            accum+=W[i]*W[j]*qi*qi*qj*qj*f[i]*f[j]*float(
                np.dot(wmu,sigma_v))
    # Two spins each, 4π×2π angles, (2π)^−6: 4/(8π⁴).
    return T**6/(2*math.pi**4)*accum


def initial_fd_density(T,m=1.,*,n=128,cut=48):
    """Number density of ONE of e− or e+ with two spins at μ=0."""
    x,w=leggauss(n)
    q=cut*(x+1)/2
    return T**3/math.pi**2*float(np.dot(cut*w/2,q*q*
                            expit(-np.sqrt(q*q+(m/T)**2))))


def quantum_event_balance(Eminus,Eplus,k1,k2,T,mumin,muplus=0.,mugamma=0.):
    """Local gain/loss on physical on-shell e−+e+<->γ+γ event.

    Spin/polarisation and actual QED phase space absent from this LOCAL
    event factor; event ledger does not give an integrated collision rate.
    """
    if T<=0 or min(Eminus,Eplus,k1,k2)<=0 or not math.isclose(
            Eminus+Eplus,k1+k2,rel_tol=1e-12):
        raise ValueError("positive energies and exact shell required")
    def f(E,mu):return float(expit(-(E-mu)/T))
    def b(E):return float(1/np.expm1((E-mugamma)/T))
    if min(k1,k2)<=mugamma: raise ValueError("boson chemical potential")
    fm,fp=f(Eminus,mumin),f(Eplus,muplus)
    b1,b2=b(k1),b(k2)
    F=fm*fp*(1+b1)*(1+b2)
    R=(1-fm)*(1-fp)*b1*b2
    affinity=(mumin+muplus-2*mugamma)/T
    return {"forward":F,"reverse":R,"ratio":F/R,
            "affinity_over_T":affinity,
            "entropy_event_factor":(F-R)*math.log(F/R)}


if __name__=="__main__":
    import argparse,json
    a=argparse.ArgumentParser()
    a.add_argument("--T-over-me",type=float,default=.05)
    a.add_argument("--n",type=int,default=32)
    o=a.parse_args()
    T=o.T_over_me
    r=initial_fd_loss_density(T,n=o.n)
    d=initial_fd_density(T)
    print(json.dumps({"T_over_me":T,"n":o.n,
      "relative_sigma_v_threshold":r/d**2/(math.pi*ALPHA**2),
      "no_photon_bose_or_inverse":True},indent=2))
