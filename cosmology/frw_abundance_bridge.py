"""Conditional relic-density normalization bridge for frozen Bose number.

A computation of n_phi/T^3 in a toy FRW model is NOT a measured relic.
Assume instantaneous entropy-decoupling at g_star_S=55; thereafter
conserved comoving number, T_gamma=2.7255 K and g_star_S(today)=3.909.
For arbitrary frozen positive PSD whose number moment is known, its
present nonrelativistic rest-mass density scales with that moment.
The shape (hence Delta Neff and P(k)) CANNOT be recovered from number.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scipy.special import zeta
from thermal_bose_class import thermal_background

N_BE = float(zeta(3.0,1.0))/math.pi**2


def relic_from_freezeout_number(n_over_T3, m_eV=86.,
                               gstarS_dec=55., gstarS_today=3.909):
    if n_over_T3<0 or m_eV<=0 or gstarS_dec<=0 or gstarS_today<=0:
        raise ValueError("Physical nonnegative abundance, positive mass/entropy required")
    anchor=thermal_background(m=m_eV,gdec=gstarS_dec,gtoday=gstarS_today)
    omega=anchor["omega_h2"]*(n_over_T3/N_BE)
    return {
        "frozen_n_phi_over_T3":float(n_over_T3),
        "mass_eV":float(m_eV),
        "gstarS_dec":float(gstarS_dec),
        "gstarS_today":float(gstarS_today),
        "thermal_Bose_number_over_T3":N_BE,
        "number_fraction_relative_thermal":float(n_over_T3/N_BE),
        "hypothetical_present_Omega_phi_h2":float(omega),
        "thermal_reference_Omega_phi_h2":float(anchor["omega_h2"]),
        "assumptions":"Decoupled stable NR relic, conserved comoving particle number, constant g* during kinetic pilot followed by separate hypothesized standard SM entropy release. No late number-changing interactions.",
        "not_inferred":"No Delta Neff, transfer P(k), Bose condensation, Lyman-alpha or likelihood prediction from the number moment alone.",
    }


def required_freezeout_number(target_omega=.12,m_eV=86.,
                              gstarS_dec=55.,gstarS_today=3.909):
    if target_omega<0:
        raise ValueError("Nonnegative target required")
    anchor=thermal_background(m=m_eV,gdec=gstarS_dec,gtoday=gstarS_today)
    return float(target_omega/anchor["omega_h2"]*N_BE)


def run(quick=False):
    kwargs=dict(m_eV=86.,gstarS_dec=55.,gstarS_today=3.909)
    target=.12
    required=required_freezeout_number(target_omega=target,**kwargs)
    baseline=relic_from_freezeout_number(N_BE,**kwargs)
    from_pilot={
        ".0001":.0000864,
        ".001":.000864,
        ".01":.008704,
        ".1":.24216,
        ".2":.36948,
    }
    # The five numbers are ROUNDED summaries of earlier independently
    # integrated nq=80 FRW toy, retained as labeled prior-result inputs;
    # re-run frw_fixed_bath_pilot.py for full-precision reproduced values.
    chosen=(".01",".1",".2") if quick else tuple(from_pilot)
    rows=[]
    for k in chosen:
        d=relic_from_freezeout_number(from_pilot[k],**kwargs)
        rows.append({"pilot_kappa_rounded":float(k),
                     "pilot_nq":80,"pilot_nphi_rounded_input":from_pilot[k],
                     "hypothetical_Omega_h2":d["hypothetical_present_Omega_phi_h2"],
                     "ratio_to_target":d["hypothetical_present_Omega_phi_h2"]/target})
    return {
        "scope":"Illustrative amplitude-only relic-density bridge, NOT cosmological validation",
        "anchor":"cosmology/thermal_bose_class.py physical constants/thermal real-scalar reference",
        "target_Omega_h2":target,
        "mass_eV":86.,
        "assumed_entropy_gstarS_dec":55.,
        "required_frozen_number_over_T3":required,
        "required_fraction_of_initial_mothers_at_x0p05":required/.3652061032773952,
        "thermal_reference":baseline,
        "rounded_prior_FRW_pilot_scenarios":rows,
        "warning":"The fixed-bath FRW pilot assumed constant g_star_S; the gS=55 -> 3.909 entropy dilution is separately HYPOTHETICAL and not dynamically computed. The mother inventory at target would suffer ~34% depletion in a closed source. No PSD-shape, Delta Neff or gravitational observables can be inferred from n alone."
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--quick",action="store_true")
    p.add_argument("--output",default="outputs/conditional_frw_relic_abundance.json")
    opts=p.parse_args()
    result=run(opts.quick)
    path=Path(opts.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
