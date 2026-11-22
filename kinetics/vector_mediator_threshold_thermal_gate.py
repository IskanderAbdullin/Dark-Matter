"""Conditional gauge-broken vector thermal and on-shell threshold audit.

Ideal massive spin-1 Bose EOS IF mediator equilibrates; no claim it
does. Independent Bessel series validates the phase-space integrals.
Decay branching widths require unknown couplings and SM mixing.
Natural units GeV. Original phi is real Z2-odd and U1-neutral.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import kv

ME_GEV = .00051099895


def mediator_ideal_equilibrium(mV_over_me: float, T_over_me: float, nodes=512):
    """Conditional massive vector equilibrium rho/P, gstar contribution.

    Requires 3 physical massive Proca polarizations and mu_V=0.
    If U1 symmetry is restored or vector never equilibrates, inapplicable.
    """
    if not (math.isfinite(mV_over_me) and math.isfinite(T_over_me)
            and mV_over_me>0 and T_over_me>0 and isinstance(nodes,int)
            and nodes>=32):
        raise ValueError("finite positive masses, T and order")
    x=mV_over_me/T_over_me
    z,w=leggauss(nodes)
    q=30*(z+1)
    w=30*w
    eps=np.hypot(q,x)
    f=1/np.expm1(eps)
    g=3
    rho=g/(2*math.pi**2)*float(np.dot(w,q*q*eps*f))
    press=g/(6*math.pi**2)*float(np.dot(w,q**4/eps*f))
    return {"mV_over_T":x,"rho_over_T4":rho,"P_over_T4":press,
            "delta_gstar_rho":30*rho/math.pi**2,
            "delta_gstar_s":45*(rho+press)/(2*math.pi**2),
            "ideal_thermal_equilibrium_assumed":True}


def independent_bessel_eos(x: float, terms=1024):
    """Independent fugacity Bessel K1/K2 massive Bose rho,P (g=3)."""
    if not math.isfinite(x) or x<=0 or terms<=0:
        raise ValueError
    rho=pres=0.
    for k in range(1,terms+1):
        z=k*x
        k1=kv(1,z);k2=kv(2,z)
        rho+=(3*k2/(k*k)+x*k1/k)
        pres+=k2/(k*k)
    pref=3*x*x/(2*math.pi**2)
    return {"rho_over_T4":pref*rho,"P_over_T4":pref*pres}



def dirac_vector_width_coefficient(mV,mchi):
    """Gamma(V->chi chibar)/(gchi*qchi)^2/mV for one Dirac fermion."""
    if mV<=0 or mchi<0:
        raise ValueError
    if mV<=2*mchi:
        return 0.
    x=mchi/mV
    return (1+2*x*x)*math.sqrt(1-4*x*x)/(12*math.pi)


def gauge_thresholds(*,m1_over_me=.1,m2_over_me=.02,
                     mphi_eV=86.,mV_over_me=.03,G_times_me2=.001):
    """Exact rest thresholds only, not thermal cross sections or widths.

    Gauge Higgs mass requires a NEW Z2-even charged complex scalar S
    with vev, OR a Stueckelberg alternative. Original real odd phi
    cannot be this symmetry-breaking charged scalar.
    """
    if not all(math.isfinite(x) and x>0 for x in
               (m1_over_me,m2_over_me,mphi_eV,mV_over_me,G_times_me2)):
        raise ValueError
    mphi_over_me=mphi_eV*1e-9/ME_GEV
    return {
      "masses_over_me":{"chi1":m1_over_me,"chi2":m2_over_me,
                         "original_real_phi":mphi_over_me,"mediator":mV_over_me},
      "mediator_keV":mV_over_me*ME_GEV*1e6,
      "contact_G_times_mV_squared":G_times_me2*mV_over_me**2,
      "yukawa_chi1_to_chi2_phi_open":m1_over_me>m2_over_me+mphi_over_me,
      "chi1_anti_chi1_to_2V_open_at_rest":m1_over_me>mV_over_me,
      "chi2_anti_chi2_to_2V_open_at_rest":m2_over_me>mV_over_me,
      "V_to_e_minus_e_plus_open_at_rest":mV_over_me>2.,
      "V_to_chi1_pair_open_at_rest":mV_over_me>2*m1_over_me,
      "V_to_chi2_pair_open_at_rest":mV_over_me>2*m2_over_me,
      "V_to_two_active_nu_flavors_open":True,
      "nu_total_tree_width_over_gprime2_mV":1/(12*math.pi),
      "chi1_pair_width_over_gchi2qchi2_mV":dirac_vector_width_coefficient(mV_over_me,m1_over_me),
      "chi2_pair_width_over_gchi2qchi2_mV":dirac_vector_width_coefficient(mV_over_me,m2_over_me),
      "mandatory_chi2_vector_current_for_real_phi_equal_charge":True,
      "not_a_thermal_reaction_rate_or_viable_model_point":True
    }


if __name__=="__main__":
    import json
    cases={}
    for x in (.03,.1,.3,1.,3.):
        cases[str(x)]={"thermo":mediator_ideal_equilibrium(x,1),
                       "thresholds":gauge_thresholds(mV_over_me=x)}
    print(json.dumps(cases,indent=2))
