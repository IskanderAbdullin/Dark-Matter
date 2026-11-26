"""Finite-charge ν+Zprime Bose-mode time evolution, restricted physical prototype.

Microscopic Zprime↔ννbar mode C_V(p) (FD Pauli, Bose inverse), integrated
as exact exponential at temporarily frozen ν FD(T,mu); reconstruct ν
T,mu each step from exact total Nν+2NV and Eν+EV conservation.
Closure ASSUMES instant number/energy-conserving ν redistribution.
It does NOT compute actual SM ν elastic/number-changing rates,
FRW, dark χ/φ kinetics or physical cosmological relic.
Dimensionless tau=Γnu,vac t; time splitting must converge.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import expit
from scipy.optimize import root
from numpy.polynomial.legendre import leggauss
from kinetics.vector_neutrino_equilibrium_charges import equilibrium_from_empty_vector


def grid(cut=35.,n=96):
    x,w=leggauss(n)
    return cut*(x+1)/2,cut*w/2


def nu_moments(T,y,*,q,w):
    f=expit(y-q)
    N=4*T**3/(2*math.pi**2)*float(np.dot(w,q*q*f))
    E=4*T**4/(2*math.pi**2)*float(np.dot(w,q**3*f))
    S=(4*E/3-y*T*N)/T
    return N,E,S


def vector_entropy(b):
    v=np.zeros_like(b)
    a=b>0
    z=b[a]
    v[a]=(1+z)*np.log1p(z)-z*np.log(z)
    return v


def coeff_thermal(m,p,T,y):
    """Exact angular Pauli ννbar coefficient for equal ν/antinu number mu."""
    E=np.hypot(p,m)
    loglo=np.logaddexp(0,y-(E-p)/(2*T))
    loghi=np.logaddexp(0,y-(E+p)/(2*T))
    a=1-2*T/p*(loglo-loghi)
    return np.where(p<1e-7*T,1-2*expit(y-m/(2*T)),a)


def integrate(*,mV_T0=.03,total_tau=1000.,steps=100,n=96,cut=35.,record=False):
    if not(mV_T0>0 and total_tau>0 and steps>=1 and n>=16 and cut>0):
        raise ValueError
    m=mV_T0
    q,w=grid(cut,n)
    p=q.copy()  # fixed PHYSICAL momentum in T0=1 units
    E=np.hypot(p,m)
    DOS=3/(2*math.pi**2)*w*p*p
    qnu,wnu=grid(50.,n)
    N0,E0,S0=nu_moments(1.,0.,q=qnu,w=wnu)
    b=np.zeros_like(p)
    T=1.;y=0.
    def state():
        Nv=float(np.dot(DOS,b));Ev=float(np.dot(DOS*E,b))
        N,E_n,S_n=nu_moments(T,y,q=qnu,w=wnu)
        Sv=float(np.dot(DOS,vector_entropy(b)))
        return {"T":T,"y":y,"Nv":Nv,"Ev":Ev,
                "N_total_relative":(N+2*Nv)/N0-1,
                "E_total_relative":(E_n+Ev)/E0-1,
                "entropy_over_initial":(S_n+Sv)/S0}
    history=[state()]
    for step in range(steps):
        a=coeff_thermal(m,p,T,y)
        if np.any(a<=0):
            raise ArithmeticError("nonpositive Pauli coefficient")
        beq=1/np.expm1(E/T-2*y)
        b=beq+(b-beq)*np.exp(-(total_tau/steps)*(m/E)*a)
        Nv=float(np.dot(DOS,b));Ev=float(np.dot(DOS*E,b))
        Ntarg=N0-2*Nv;Etarg=E0-Ev
        if Ntarg<=0 or Etarg<=0:
            raise ArithmeticError("step too large or lost bath")
        def residual(z):
            Tnew=math.exp(z[0]);yn=z[1]
            Nnew,En,_=nu_moments(Tnew,yn,q=qnu,w=wnu)
            return (Nnew/Ntarg-1,En/Etarg-1)
        sol=root(residual,[math.log(T),y],tol=1e-11)
        if not sol.success and max(abs(x) for x in residual(sol.x))>1e-10:
            raise RuntimeError(sol.message)
        T=math.exp(sol.x[0]);y=float(sol.x[1])
        current=state()
        if current["entropy_over_initial"]+1e-11<history[-1]["entropy_over_initial"]:
            raise ArithmeticError("entropy decrease; reduce dtau")
        history.append(current)
    eq=equilibrium_from_empty_vector(m,n=256)
    final=history[-1]
    return {"mV_T0":m,"total_tau_Gamma0_t":total_tau,"substeps":steps,
            "vector_modes":n,"momentum_cut_T0":cut,
            "initial":history[0],"final":final,
            "reference_closed_equilibrium":eq,
            "delta_T_to_equilibrium":final["T"]-eq["T_final_over_T0"],
            "delta_y_to_equilibrium":final["y"]-eq["mu_nu_over_T_final"],
            "entropy_monotonic":True,
            "max_abs_N_relative_error":max(abs(h["N_total_relative"]) for h in history),
            "max_abs_E_relative_error":max(abs(h["E_total_relative"]) for h in history),
            "requires_instantaneous_nu_thermal_number_conserving_closure":True,
            "no_real_SM_weak_or_FRW_or_chi":True,
            "history":history if record else None}


if __name__=="__main__":
    import json,argparse
    p=argparse.ArgumentParser()
    p.add_argument("--steps",type=int,default=100)
    p.add_argument("--tau",type=float,default=1000.)
    p.add_argument("--modes",type=int,default=96)
    o=p.parse_args()
    print(json.dumps(integrate(total_tau=o.tau,steps=o.steps,n=o.modes),indent=2))
