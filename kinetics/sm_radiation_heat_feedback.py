"""Finite-heat-capacity Standard Model RADIATION bath coupled to elastic chi1 BGK.

This is the next honest conservation stage beyond silently resetting
chi1 to FD(T). One fixed-a collision substep:
  1. number-project the mother's target FD onto the actual Nchi1;
  2. take its exactly positive Pauli-safe exponential BGK step;
  3. transfer the *actual* dark comoving energy gain/loss to/from SM;
  4. derive the NEW bath temperature from finite radiation heat capacity:
       a^4 rho_SM = pi²/30 * gstar * (a*T_SM)^4.
This treats SM as internally thermally equilibrated massless radiation
with one constant gstar_s=gstar_rho, not a microscopic 2->2 scattering
model and not realistic changing gstar(T), finite-mass threshold or
self-consistent Friedmann evolution.

Entropy test: S_com,SM=a³ * 2pi²/45*gstar*T³,
              S_com,chi=-g/(2pi²)int dq q² [f ln f+(1-f)ln(1-f)].
If the finite-step entropy falls below tolerance, fail closed: the
held-T bath approximation was too coarse and MUST be subdivided.

This module is a conservatively coupled collision *substep* at fixed a.
Full FRW evolution still requires physical pressure work separately.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import xlogy
from number_conserving_sm_bath import (
    exact_frozen_bath_step,_weights_and_state
)


def comoving_radiation_energy(a,T,gstar=55.):
    if (not all(math.isfinite(v) for v in (a,T,gstar)) or
            a<=0 or T<=0 or gstar<=0):
        raise ValueError("Positive finite scale factor, bath T and radiation dof required")
    return float(math.pi**2/30*gstar*(a*T)**4)


def comoving_radiation_entropy(a,T,gstar=55.):
    _=comoving_radiation_energy(a,T,gstar)
    return float(2*math.pi**2/45*gstar*(a*T)**3)


def _comoving_fermion_entropy(q,f,g):
    _,occupation,W=_weights_and_state(q,f,g)
    x=np.clip(occupation,1e-15,1-1e-15)
    return float(np.dot(W,-xlogy(x,x)-xlogy(1-x,1-x)))


def radiation_elastic_exchange_step(q,f,a=1.18,mass=.1,
                                     T_bath=.82,g1=4.,gamma_el=1.2,
                                     dt=.4,gstar=55.):
    """Exact conservative energy transfer with a finite radiation bath.

    Does NOT assert that mu=0, the bath is infinite, nor that Γ_el
    comes from a particular Standard Model scattering matrix element.
    """
    E_bath_before=comoving_radiation_energy(a,T_bath,gstar)
    S_bath_before=comoving_radiation_entropy(a,T_bath,gstar)
    S_mother_before=_comoving_fermion_entropy(q,f,g1)
    molecular=exact_frozen_bath_step(
        q,f,a=a,mass=mass,Tbath=T_bath,
        g=g1,gamma=gamma_el,dt=dt)
    qarr=np.asarray(q,dtype=float)
    before=np.asarray(f,dtype=float)
    norm=g1/(2*math.pi**2)
    Ecom=np.hypot(qarr,a*mass)
    E_mother_before=float(norm*np.trapezoid(qarr**2*Ecom*before,qarr))
    E_mother_after=float(norm*np.trapezoid(
        qarr**2*Ecom*molecular["after"],qarr))
    E_dark_delta=E_mother_after-E_mother_before
    E_bath_after=E_bath_before-E_dark_delta
    if not math.isfinite(E_bath_after) or E_bath_after<=0:
        raise ValueError("Inadequate radiation bath heat capacity for energy exchange")
    T_after=float((E_bath_after/
                   (math.pi**2/30*gstar*a**4))**.25)
    S_after=(comoving_radiation_entropy(a,T_after,gstar)+
             _comoving_fermion_entropy(q,molecular["after"],g1))
    S_before=S_bath_before+S_mother_before
    entropy_delta=S_after-S_before
    if entropy_delta<-1e-9*max(abs(S_before),1.):
        raise ArithmeticError("Finite-step BGK held-bath approximation reduces total entropy: subdivide")
    result={
        "q":qarr,
        "mother_before":before,"mother_after":molecular["after"],
        "equilibrium_target_at_T_before":molecular["target_fermi_dirac"],
        "mother_number_before":molecular["number_before"],
        "mother_number_after":molecular["number_after"],
        "mother_number_change":molecular["comoving_number_residual"],
        "mother_projected_mu":molecular["chemical_mu_number_projected"],
        "mother_energy_before":E_mother_before,
        "mother_energy_after":E_mother_after,
        "SM_radiation_comoving_energy_before":E_bath_before,
        "SM_radiation_comoving_energy_after":E_bath_after,
        "total_energy_before":E_mother_before+E_bath_before,
        "total_energy_after":E_mother_after+E_bath_after,
        "total_comoving_collision_energy_residual":
            E_mother_after+E_bath_after-E_mother_before-E_bath_before,
        "T_SM_before":T_bath,
        "T_SM_after":T_after,
        "gstarS_and_gstarRho_both_constant_and_equal":gstar,
        "mother_entropy_before":S_mother_before,
        "SM_entropy_before":S_bath_before,
        "total_entropy_change":entropy_delta,
        "fermionic_KL_before":molecular["fermionic_KL_before"],
        "fermionic_KL_after":molecular["fermionic_KL_after"],
        "Pauli_occupation_min_after":molecular["minimum_occupation"],
        "Pauli_occupation_max_after":molecular["maximum_occupation"],
        "finite_temperature_feedback_implemented":True,
        "time_integrator_first_order_bath_feedback":True,
        "is_one_fixed_a_collision_substep":True,
        "microphysical_SM_scattering_matrix_element_derived":False,
        "full_FRW_pressure_and_1to2_scatterings_integrated":False,
        "physical_86eV_relic_or_CLASS_ready":False
    }
    return result


def demonstration():
    from scipy.special import expit
    q=np.geomspace(1e-4,30,96)
    initial=np.clip(expit(-np.hypot(q,.1))*
                    (1+.05*np.tanh(np.log1p(q))),0,1)
    r=radiation_elastic_exchange_step(q,initial)
    return {k:v for k,v in r.items() if not isinstance(v,np.ndarray)}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="outputs/finite_SM_radiation_thermal_contact.json")
    opt=parser.parse_args()
    r=demonstration()
    path=Path(opt.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(r,indent=2),encoding="utf-8")
    print(json.dumps(r,indent=2))
