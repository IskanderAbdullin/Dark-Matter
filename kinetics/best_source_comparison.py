"""Fail-closed independent BEST vacuum-source Monte Carlo comparison.

INPUT: one or more JSON files actually emitted by
best_external_vacuum_reference.py --execute, which checks Git SHA.
This validator NEVER calls BEST, and cannot authenticate a fabricated
JSON file. It checks analytic rate, three species stoichiometry,
one or more grids and stochastic sampling settings; with insufficient
independent settings it labels the result UNVERIFIED.

Even a passed instantaneous vacuum-rate comparison does NOT validate
Bose-stimulated nonlinear production, zero-mode, FRW, or CLASS.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from best_external_vacuum_reference import BEST_PIN, benchmark_protocol


REQUIRED_SPECIES=("chi1","phi","chi2")


def inspect_record(data,rate_tolerance=.1,charge_tolerance=.1):
    if not isinstance(data,dict) or not data.get("external_solver_executed"):
        raise ValueError("Unexecuted or invalid external BEST reference record")
    if data.get("BEST_source_commit_verified")!=BEST_PIN or not data.get(
            "BEST_git_revision_matched_pinned_protocol"):
        raise ValueError("Externally executed BEST record does not match pinned Git SHA")
    proto=data.get("model_and_protocol",{})
    required=benchmark_protocol(
        mass_over_T=proto.get("parent_mass_over_T",.1),
        parent_dof=proto.get("parent_dof",4),
        rest_width=proto.get("rest_width",1.),
        grid_n=proto.get("grid_n",100),
        qmin=proto.get("qmin",1e-4),qmax=proto.get("qmax",30.),
        vegas_neval=proto.get("vegas_neval_per_grid_point",20000))
    if (proto.get("pinned_best_git_commit")!=BEST_PIN or
            proto.get("daughter_bose_dof")!=1 or
            proto.get("daughter_fermi_dof")!=4 or
            abs(proto.get("best_spin_summed_matrix_element_squared",0)-
                required["best_spin_summed_matrix_element_squared"])>
                1e-10*required["best_spin_summed_matrix_element_squared"]):
        raise ValueError("Mismatched independent physical matrix-element conventions")
    numeric=data.get("instantaneous_integrated_number_rates",{})
    if not all(k in numeric for k in REQUIRED_SPECIES):
        raise ValueError("Missing one or more independent species rate moments")
    rates={k:float(numeric[k]) for k in REQUIRED_SPECIES}
    eraw=data.get("instantaneous_integrated_energy_rates",{})
    if not all(k in eraw for k in REQUIRED_SPECIES):
        raise ValueError("Missing one or more independent species ENERGY rate moments")
    energy_rates={k:float(eraw[k]) for k in REQUIRED_SPECIES}
    if not all(math.isfinite(v) for v in energy_rates.values()):
        raise ValueError("Non-finite external ENERGY rate data")
    if not all(math.isfinite(v) for v in rates.values()):
        raise ValueError("Non-finite external rate data")
    target=required["independently_derived_vacuum_source_n_over_T3_per_time"]
    if rate_tolerance<=0 or charge_tolerance<=0:
        raise ValueError("Positive diagnostic tolerances required")
    phi=abs(rates["phi"]);chi2=abs(rates["chi2"]);chi1=abs(rates["chi1"])
    scale=max(phi,chi2,chi1,1e-300)
    source_error=abs(rates["phi"]/target-1)
    stoichiometry=max(abs(rates["chi1"]+rates["phi"]),
                     abs(rates["chi1"]+rates["chi2"]),
                     abs(rates["phi"]-rates["chi2"]))/scale
    signs=rates["chi1"]<0 and rates["phi"]>0 and rates["chi2"]>0
    e_target=required["independently_derived_total_daughter_energy_gain_over_T4_per_time"]
    e_scale=max(abs(x) for x in energy_rates.values())
    e_scale=max(e_scale,1e-300)
    energy_balance=abs(sum(energy_rates.values()))/e_scale
    energy_reference_error=abs(-energy_rates["chi1"]/e_target-1)
    energy_signs=(energy_rates["chi1"]<0 and energy_rates["phi"]>0 and
                  energy_rates["chi2"]>0)
    energy_pass=energy_signs and energy_balance<=charge_tolerance and (
        energy_reference_error<=rate_tolerance)
    return {
        "BEST_git_commit_verified":BEST_PIN,
        "mass_over_T":required["parent_mass_over_T"],
        "rest_width":required["rest_width"],
        "parent_dof":required["parent_dof"],
        "qmin":required["qmin"],
        "qmax":required["qmax"],
        "grid_n":required["grid_n"],
        "vegas_neval":required["vegas_neval_per_grid_point"],
        "raw_rates":rates,
        "raw_physical_energy_rates":energy_rates,
        "analytic_vacuum_source":target,
        "independent_exact_parent_energy_loss":e_target,
        "relative_parent_energy_loss_error":-energy_rates["chi1"]/e_target-1,
        "relative_total_energy_rate_residual":energy_balance,
        "energy_rate_signs_correct":energy_signs,
        "physical_energy_rate_checks_pass":energy_pass,
        "relative_phi_rate_error":rates["phi"]/target-1,
        "relative_three_charge_stoichiometry_error":stoichiometry,
        "correct_physical_signs":signs,
        "initial_source_within_rate_tolerance":source_error<=rate_tolerance,
        "all_three_numbers_within_charge_tolerance":
            stoichiometry<=charge_tolerance,
        "single_record_passes_preliminary_criteria":
            signs and source_error<=rate_tolerance and stoichiometry<=charge_tolerance and energy_pass,
        "warning":"A single stochastic BEST rate pass does not imply independent convergence or physically correct late-time dark matter."
    }


def evaluate_records(records,rate_tolerance=.1,charge_tolerance=.1):
    if not records:
        return {
            "status":"UNVERIFIED: no independently executed BEST results provided",
            "records_checked":0,
            "source_agreement_verified":False,
            "freezeout_or_CLASS_validated":False,
            "checks":[]
        }
    checks=[inspect_record(r,rate_tolerance,charge_tolerance) for r in records]
    grid_values=sorted({int(r["grid_n"]) for r in checks})
    neval_values=sorted({int(r["vegas_neval"]) for r in checks})
    at_least_two_grids=len(grid_values)>=2
    at_least_two_neval=len(neval_values)>=2
    all_success=all(x["single_record_passes_preliminary_criteria"] for x in checks)
    same_physics=all((x["mass_over_T"],x["rest_width"],x["parent_dof"],
                      x["qmin"],x["qmax"])==(
                      checks[0]["mass_over_T"],checks[0]["rest_width"],
                      checks[0]["parent_dof"],checks[0]["qmin"],
                      checks[0]["qmax"]) for x in checks)
    preliminary=all_success and same_physics and at_least_two_grids and (
        at_least_two_neval) and len(checks)>=4
    return {
       "status":("PRELIMINARY SOURCE RATE AGREEMENT OF SUPPLIED BEST RECORDS"
                 if preliminary else "UNVERIFIED OR FAILING INDEPENDENT RATE EVIDENCE"),
       "records_checked":len(checks),
       "q_grid_sizes":grid_values,"vegas_samples_per_grid_values":neval_values,
       "same_physical_decay_scenario_across_all_runs":same_physics,
       "multiple_independent_momentum_grids":at_least_two_grids,
       "multiple_monte_carlo_sample_counts":at_least_two_neval,
       "all_pass_individual_rate_and_number_checks":all_success,
       "source_agreement_verified":preliminary,
       "freezeout_or_CLASS_validated":False,
       "checks":checks,
       "warning":"Files are externally supplied data; content/provenance cannot be cryptographically authenticated by this comparator. BEST runner verifies source SHA. A pass assesses only initial thermal vacuum rate, not MC uncertainties, nonlinear Bose, zero-mode or cosmology."
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("inputs",nargs="*",help="BEST JSON outputs from separate --execute runs")
    p.add_argument("--output",default="outputs/best_source_comparison.json")
    p.add_argument("--rate-tolerance",type=float,default=.1)
    p.add_argument("--charge-tolerance",type=float,default=.1)
    args=p.parse_args()
    records=[json.loads(Path(file).read_text(encoding="utf-8"))
             for file in args.inputs]
    report=evaluate_records(records,args.rate_tolerance,args.charge_tolerance)
    dest=Path(args.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
