"""Independent convergence audit: static mother moments != nonlinear source moments.

The THREE physical collision invariants are verified for every bridged
on-shell mesh. Run a short fixed-T, dynamically CLOSED 1<->2 reaction:
  Nq=199...299, Nparent=5099...10199, M1/T=0.1, Gamma*t=0.001.
Compare upper daughter cutoff and time-integration tolerances. This
does NOT establish primordial freeze-out, the full shape of f(q),
a zero-mode/condensate, or cosmological dark-matter viability.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from schur_closed_bdf import solve_closed_schur
from sparse_kinetic_jacobian import SparseClosedDecay
from resonant_graph_repair import (
    bridged_indices,exact_mother_thermal_moments,parent_grid_moments
)


def one_run(max_index, rtol=1e-5, atol=1e-11, tend=.001):
    mesh=bridged_indices(max_index=max_index)
    model=SparseClosedDecay(mesh)
    assert model.graph["active_graph_collision_invariants"]==3
    before=parent_grid_moments(mesh)
    reference=exact_mother_thermal_moments()
    errors={key:float(before[key]/reference[key]-1)
            for key in ("number_over_T3","energy_over_T4","pressure_over_T4")}
    result=solve_closed_schur(model,tend=tend,rtol=rtol,atol=atol)
    return {
      "max_index":max_index,
      "daughter_modes":model.nq,
      "mother_modes":model.nparent,
      "rtol":rtol,"atol":atol,"t_Gamma":tend,
      "initial_parent_FD_relative_errors":errors,
      "dynamic_phi_number_over_T3":result["n_phi_over_T3"],
      "dynamic_phi_mean_q_over_T":result["mean_phi_q"],
      "dynamic_phi_fraction_below_q0p1":result["phi_fraction_below_q0p1"],
      "max_bose_occupation_grid_dependent":
          result["max_boson_occupation_grid_dependent"],
      "maximum_absolute_charge_or_energy_drift":max(
          abs(v) for v in result["number_and_energy_residuals"].values()),
      "integration_nfev":result["nfev"],
      "graph_invariants":result["graph_collision_invariants"],
    }


def relative_delta(current,base):
    return float(current/base-1)


def run(quick=False):
    cutoffs=(5000,7000) if quick else (5000,7000,8500,10000)
    results=[one_run(m) for m in cutoffs]
    tighter_index=7000 if quick else 8500
    tightened=one_run(tighter_index,rtol=1e-6,atol=1e-12)
    original=next(r for r in results if r["max_index"]==tighter_index)
    time_number_difference=abs(relative_delta(
        original["dynamic_phi_number_over_T3"],
        tightened["dynamic_phi_number_over_T3"]))
    late_a,late_b=results[-2:]
    cutoff_number_difference=abs(relative_delta(
        late_b["dynamic_phi_number_over_T3"],
        late_a["dynamic_phi_number_over_T3"]))
    late_momentum_difference=abs(relative_delta(
        late_b["dynamic_phi_mean_q_over_T"],
        late_a["dynamic_phi_mean_q_over_T"]))
    return {
       "title":"Cutoff-vs-time error audit of invariant-clean short-time closed quantum kinetics",
       "scope":"M1/T=.1, fixed-T flat space, t*Gamma=.001, initial FD mother, empty daughters",
       "status":"PRELIMINARY, not a cosmologically frozen spectrum",
       "cutoff_variations":results,
       "time_refinement_at_max_index":tighter_index,
       "time_refinement_result":tightened,
       "relative_number_time_discrepancy":time_number_difference,
       "relative_number_cutoff_discrepancy_last_two":cutoff_number_difference,
       "relative_mean_q_cutoff_discrepancy_last_two":late_momentum_difference,
       "ratio_cutoff_discrepancy_to_time_discrepancy":
           cutoff_number_difference/max(time_number_difference,1e-15),
       "initial_FD_moments_are_not_proxy_for_dynamic_source_accuracy":True,
       "warnings":"No independent convergence in delta, fine-core size, stride, lower IR cutoff or final time. No zero-mode physics, FRW, physical freeze-out or CLASS. The 1<->2 collision kernel is still finite-cell approximate.",
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick",action="store_true")
    parser.add_argument("--output",default="outputs/closed_cutoff_and_time_audit.json")
    arg=parser.parse_args()
    result=run(arg.quick)
    p=Path(arg.output);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
