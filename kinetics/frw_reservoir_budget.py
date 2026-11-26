"""Comoving mother-reservoir budget diagnostic for fixed-thermal-bath FRW pilot.

This is NOT a replacement kinetic calculation. Given a freezeout result from
frw_fixed_bath_pilot, test the exact upper bound from the independently
conserved closed reaction charge:
 a^3(n_chi1+n_phi) = const.
For initially empty phi and chi2 and g*S constant, any CLOSED 1<->2
evolution obeys n_phi/T^3 <= n_chi1(initial)/T_initial^3.

If a fixed-bath pilot breaches the bound, external replenishment is
mathematically mandatory to interpret that pilot. If it stays below,
large implied depletion can still invalidate fixed-mother dynamics.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit, zeta

from frw_fixed_bath_pilot import FRWFixedBathPilot


def thermal_mother_number(x_initial: float, g1: int = 4) -> float:
    """Dimensionless thermal Fermi density n/T^3 for mass M/T=x_initial."""
    if x_initial < 0 or g1 <= 0:
        raise ValueError("physical x>=0 and g1>0 required")
    answer = quad(lambda q:q*q*expit(-math.hypot(q,x_initial)),
                  0,np.inf,epsabs=1e-12,epsrel=1e-12,limit=200)[0]
    return g1/(2*math.pi**2)*answer


def mother_budget(nphi_over_T3: float, x_initial=.05, g1=4) -> dict:
    if nphi_over_T3 < 0:
        raise ValueError("Negative produced density")
    parent = thermal_mother_number(x_initial,g1)
    massless_parent = (3./4.)*g1*float(zeta(3,1))/math.pi**2
    consumed = nphi_over_T3/parent
    return {
        "initial_x_M1_over_T":x_initial,
        "initial_thermal_parent_n_over_T3":parent,
        "massless_parent_upper_limit_n_over_T3":massless_parent,
        "produced_phi_n_over_T3":nphi_over_T3,
        "closed_reaction_minimum_parent_depletion_fraction":consumed,
        "closed_without_external_parent_source_possible_by_number":consumed<=1.,
        "minimum_required_external_net_parent_replenishment_over_T3":
            max(0.,nphi_over_T3-parent),
        "note":"Parent-number-only necessary condition, NOT a sufficient test of a closed model; external SM processes could supply parents but must then be explicitly budgeted.",
    }


def run(quick=False):
    kappas=[.01,.2] if quick else [.0001,.001,.01,.05,.1,.2]
    cases=[]
    for kappa in kappas:
        model=FRWFixedBathPilot(nq=80,kappa=kappa)
        hist=model.integrate()
        budget=mother_budget(hist["final"]["nphi_over_T3"],hist["integration"]["x_start"])
        cases.append({"kappa":kappa,"nq":80,"freezeout":hist["final"],
                      "number_budget":budget,
                      "late_number_drift":hist["relative_number_change_from_mid_to_end"]})
    return {
        "description":"Audit of necessary parent-number budget for a hypothetical CLOSED interpretation of the FIXED-bath FRW pilot",
        "closed_reaction_charge":"a^3(n_chi1+n_phi)=constant when phi initially empty and no other reactions",
        "initial_condition":"thermal mother at x=.05, empty daughters, constant g_star_s",
        "warning":"Fixed-bath pilot explicitly allows unlimited external replenishment and does not violate its OWN assumptions. Bound exposes where a closed reinterpretation is impossible, and how much mother would need to be depleted even if bound holds.",
        "cases":cases,
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick",action="store_true")
    parser.add_argument("--output",default="outputs/frw_reservoir_budget.json")
    options=parser.parse_args()
    answer=run(quick=options.quick)
    path=Path(options.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(answer,indent=2),encoding="utf-8")
    print(json.dumps(answer,indent=2))
