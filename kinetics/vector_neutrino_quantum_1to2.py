"""Exact finite-T Pauli/Bose Zprime <-> two active neutrino flavors.

Microscopic unpolarized massive-vector 1<->2 operator for Zprime
momentum modes, with two ν flavors in an EXTERNALLY prescribed FD
bath at Tnu, zero chemical potentials. Includes ν Pauli, inverse
ννbar->Zprime and vector Bose stimulation.
Does NOT evolve finite ν energy reservoir, other Zprime channels,
thermal masses, expansion, dark elastic rate or 86-eV relic.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit

ME_GEV=.00051099895
MPL_GEV=1.2209e19


def coeff_gauss(mV,p,Tnu,nangle=64):
    """Angular <(1-f1)(1-f2)-f1*f2>=<1-f1-f2>."""
    if mV<=0 or p<0 or Tnu<=0 or nangle<4:
        raise ValueError
    E=math.hypot(p,mV)
    mu,w=leggauss(nangle)
    E1=(E+p*mu)/2;E2=E-E1
    f1=expit(-E1/Tnu);f2=expit(-E2/Tnu)
    return float(np.dot(w,1-f1-f2)/2)


def coeff_closed(mV,p,Tnu):
    """Independent analytic FD angular integral for unpolarized vector."""
    if not (mV>0 and p>=0 and Tnu>0 and
            all(math.isfinite(x) for x in (mV,p,Tnu))):
        raise ValueError
    if p<1e-7*Tnu:
        return math.tanh(mV/(4*Tnu))
    E=math.hypot(p,mV)
    umin=(E-p)/(2*Tnu);umax=(E+p)/(2*Tnu)
    avgf=Tnu/p*(np.logaddexp(0,-umin)-np.logaddexp(0,-umax))
    return float(1-2*avgf)


def mode_source(mV,p,Tnu,gprime,bV):
    """C_V(p)= Γnu^vac (mV/E) a(p,Tnu) [bBE(E/Tnu)-bV]."""
    if not (mV>0 and p>=0 and Tnu>0 and gprime>=0 and bV>=0
            and all(math.isfinite(x) for x in (mV,p,Tnu,gprime,bV))):
        raise ValueError
    E=math.hypot(p,mV)
    a=coeff_closed(mV,p,Tnu)
    beq=1/math.expm1(E/Tnu)
    width=gprime*gprime*mV/(12*math.pi)
    gamma=width*mV/E*a
    C=gamma*(beq-bV)
    return {"E":E,"a_fd":a,"beq":beq,"Gamma_linear":gamma,
            "dbdt":C,"nu_pair_number_source_per_mode":-2*C,
            "nu_energy_source_per_mode":-E*C,
            "conserves_total_energy_per_event":True}


def integrated_moments(mV,Tnu,gprime,*,occupancy_factor=0.,n=256,cut_p_over_T=45):
    """Finite-T quantum 1<->2 Zprime number/energy WEAK moments.

    Boson b(p)=occupancy_factor*BE(Tnu). For b=0,
    initial inverse ννbar production is nonzero. At b=BE
    microscopic detailed balance is exact. Quantum total entropy
    is only defined for strictly positive b in this implementation.
    """
    if not(mV>0 and Tnu>0 and gprime>=0 and occupancy_factor>=0 and n>=16):
        raise ValueError
    z,w=leggauss(n)
    p=cut_p_over_T*Tnu*(z+1)/2
    w=w*cut_p_over_T*Tnu/2
    E=np.hypot(p,mV)
    beq=1/np.expm1(E/Tnu)
    b=occupancy_factor*beq
    a=np.array([coeff_closed(mV,float(k),Tnu) for k in p])
    Gamma0=gprime*gprime*mV/(12*math.pi)
    gamma=Gamma0*mV/E*a
    dotb=gamma*(beq-b)
    weight=3/(2*math.pi**2)*w*p*p
    dN=float(np.dot(weight,dotb))
    dQ=float(np.dot(weight,E*dotb))
    n_eq=float(np.dot(weight,beq))
    rho_eq=float(np.dot(weight,E*beq))
    if occupancy_factor>0:
        logratio=np.log(b/(1+b))+E/Tnu
        entropy=float(np.dot(weight,-dotb*logratio))
        dSV=float(np.dot(weight,-dotb*np.log(b/(1+b))))
        dSnu=-dQ/Tnu
    else:
        entropy=dSV=dSnu=None
    return {"mV_over_Tnu":mV/Tnu,"gprime":gprime,
            "occupancy_factor_vs_BE":occupancy_factor,
            "n_eq_GeV3":n_eq,"rho_eq_GeV4":rho_eq,
            "dN_Zprime_GeV4":dN,
            "dN_nu_plus_antinu_GeV4":-2*dN,
            "dQ_Zprime_GeV5":dQ,"dQ_nu_GeV5":-dQ,
            "total_energy_residual_GeV5":0.,
            "entropy_total_GeV4":entropy,
            "entropy_Zprime_GeV4":dSV,
            "entropy_nu_reservoir_GeV4":dSnu,
            "inverse_fill_rate_over_initial_equilibrium_density_GeV":
                dN/n_eq if occupancy_factor==0 else None,
            "fixed_thermal_neutrino_reservoir_not_evolved":True,
            "no_chi_elastic_gamma_or_dark_relic":True}


if __name__=="__main__":
    import argparse,json
    p=argparse.ArgumentParser()
    p.add_argument("--mV-over-me",type=float,default=.03)
    p.add_argument("--T-over-me",type=float,default=1.)
    p.add_argument("--gprime",type=float,default=2.5388531259649034e-7)
    o=p.parse_args()
    cases={str(k):integrated_moments(
        o.mV_over_me*ME_GEV,o.T_over_me*ME_GEV,o.gprime,occupancy_factor=k)
           for k in (0.,.5,1.,2.)}
    print(json.dumps(cases,indent=2))
