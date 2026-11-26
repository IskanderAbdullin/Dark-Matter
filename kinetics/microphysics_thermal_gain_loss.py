"""Tree-level thermal quantum factors and physical -> BEST amplitude bridge.

This verifies on-shell BE/FD detailed balance from a concrete scalar Yukawa
vertex. It is NOT the full finite-T collision integral Gamma(T), nor a
plasma thermalisation/zero-mode calculation.
"""
from __future__ import annotations
import math
from microphysics_portal import (scalar_yukawa_width,best_spin_summed_m2)


def _positive_finite(name,x):
    if not math.isfinite(x) or x<=0:
        raise ValueError(name+" must be finite and positive")


def fermi(E,T,mu=0.):
    _positive_finite("T",T)
    if not math.isfinite(E) or E<0 or not math.isfinite(mu):
        raise ValueError("Finite E>=0 and finite chemical potential required")
    x=(E-mu)/T
    if x>=0:
        em=math.exp(-x)
        return em/(1+em)
    return 1/(1+math.exp(x))


def bose(E,T,mu=0.):
    _positive_finite("T",T)
    if not math.isfinite(E) or E<0 or not math.isfinite(mu):
        raise ValueError("Finite E>=0 and finite chemical potential required")
    x=(E-mu)/T
    if x<=0:
        raise ValueError("Bose chemical potential must be strictly below E")
    if x>700:
        return math.exp(-x)
    return 1/math.expm1(x)


def gain_loss_on_shell(E1,E2,Ephi,T,mu1=0.,mu2=0.,muphi=0.):
    """Quantum 1 <-> 2 occupation factors, no angular/phase-space integral.

    The real self-conjugate phi's chemical potential is mu_phi=0 for
    physical number-nonconserving equilibrium; nonzero mu_phi is a
    diagnostic for a formal momentarily constrained distribution.
    """
    _positive_finite("T",T)
    if (not math.isfinite(E1) or not math.isfinite(E2) or
            not math.isfinite(Ephi) or min(E1,E2,Ephi)<0 or
            not math.isclose(E1,E2+Ephi,rel_tol=1e-12,abs_tol=1e-12)):
        raise ValueError("Physical on-shell E1=E2+Ephi required")
    f1=fermi(E1,T,mu1)
    f2=fermi(E2,T,mu2)
    b=bose(Ephi,T,muphi)
    decay=f1*(1-f2)*(1+b)
    inverse=(1-f1)*f2*b
    affinity=(mu1-mu2-muphi)/T
    # Important stable near equilibrium: gain*(1-exp(-affinity)).
    if affinity < -690:
        raise ValueError("Extreme chemical affinity outside stable range")
    net=decay-inverse
    stable=decay*(-math.expm1(-affinity))
    return {"decay_gain_factor":decay,"inverse_gain_factor":inverse,
            "net_factor":net,"chemical_affinity_over_T":affinity,
            "detailed_balance_stable_net":stable,
            "occupations":{"f1":f1,"f2":f2,"bphi":b},
            "is_full_finite_temperature_rate":False}


def yukawa_best_bridge(m1_GeV,m2_GeV,mphi_GeV,T_GeV,y,parent_dof=4):
    """Dimensionless BEST units for the *same* fixed-coupling Yukawa theory.

    BEST baseline works with masses M/T and evolution time tau=t*T.
    For dimensionful |M|^2 [GeV^2], pass |M|^2/T^2 to those units.
    Includes both chi1 charge states for g1=4. NOT an executed BEST run.
    """
    _positive_finite("T",T_GeV)
    width=scalar_yukawa_width(m1_GeV,m2_GeV,mphi_GeV,y)
    if not width>0:
        raise ValueError("Open decay with nonzero Yukawa coupling required")
    amp=best_spin_summed_m2(m1_GeV,m2_GeV,mphi_GeV,y,parent_dof)
    return {
      "model":"real scalar Yukawa, Z2 chi1/phi odd, tree vacuum rate",
      "external_BEST_executed":False,
      "physical_T_GeV":T_GeV,
      "physical_Gamma_rest_GeV":width,
      "BEST_Gamma_over_T":width/T_GeV,
      "BEST_m1_over_T":m1_GeV/T_GeV,
      "BEST_m2_over_T":m2_GeV/T_GeV,
      "BEST_mphi_over_T":mphi_GeV/T_GeV,
      "BEST_spin_summed_m2_over_T2":amp/T_GeV**2,
      "BEST_parent_dof":parent_dof,
      "temperature_dependent_Gamma_derived":False,
      "microphysical_elastic_BGK_calibrated":False,
      "CLASS_ready":False,
    }
