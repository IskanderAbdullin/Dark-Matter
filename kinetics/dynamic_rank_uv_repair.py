"""Repair a changing on-shell 1<->2 graph across FRW by expanding UV support.

For each instantaneous mass ratio x=M1/T, the same sparse daughter's
fixed-q comb may lose reaction connectivity and acquire fake conserved
charges. Search increasing upper cutoffs, and accept ONLY a graph with
exactly 3 collision invariants and no isolated modes. Then INDEPENDENTLY
test its truncated thermal mother number/energy/pressure against the
continuum FD reference. Having a clean graph is necessary, NOT enough
for correct collision physics or a CLASS primordial f(q).

This is a mesh DESIGN / feasibility scan, not an evolving FRW solver.
Extending the daughter cutoff over time also changes available phase
space; physical population remapping must conserve number and energy
before any actual collision+expansion integration can use it.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from collision_invariant_rank import audit_integer_grid
from resonant_graph_repair import (
    bridged_indices, exact_mother_thermal_moments, parent_grid_moments
)


DEFAULT_UPPER=(300,400,500,700,1000,1500,2000,2500,3000,3500,4000,5000)


def repaired_instantaneous_grid(x,upper_candidates=DEFAULT_UPPER,
                                core=100,stride=50,delta=.001):
    if x<=0 or not upper_candidates:
        raise ValueError("Positive mass and at least one UV cutoff required")
    attempts=[]
    for maximum in upper_candidates:
        indices=bridged_indices(core=core,stride=stride,max_index=int(maximum))
        graph=audit_integer_grid(indices,delta=delta,mass=x)
        clean=bool(graph["accept_for_unique_chemical_equilibrium"])
        item={"upper_index":int(maximum),
              "daughter_bins":len(indices),
              "collision_invariants":graph["active_graph_collision_invariants"],
              "extra_invariants":graph["additional_invariants_beyond_number_number_energy"],
              "isolated_daughters":graph["isolated_boson_modes"]+graph["isolated_chi2_modes"],
              "graph_pass":clean}
        attempts.append(item)
        if clean:
            mom=parent_grid_moments(indices,delta=delta,mass=x)
            exact=exact_mother_thermal_moments(mass=x)
            errors={key:float(mom[key]/exact[key]-1)
                    for key in ("number_over_T3","energy_over_T4","pressure_over_T4")}
            return {
                "instantaneous_M1_over_T":x,
                "selected_upper_index":int(maximum),
                "selected_daughter_modes":len(indices),
                "selected_parent_modes":mom["active_mother_modes"],
                "selected_daughter_max_energy_over_T":
                    (max(indices)+.5)*delta,
                "number_of_collision_invariants":graph["active_graph_collision_invariants"],
                "attempts":attempts,
                "mother_thermal_relative_errors":errors,
                "independent_parent_moment_accuracy_better_than_one_percent":
                    max(abs(y) for y in errors.values())<.01,
                "collision_graph_pass":True,
                "frozen_DM_spectrum_available":False,
                "warning":"Graph connectivity repaired by UV extension, but source integration, physical remapping between grids and nonlinear/FRW convergence NOT verified."
            }
    raise ValueError(
        f"No invariant-clean on-shell grid at M/T={x} within {tuple(upper_candidates)}; "
        "do not continue FRW with fake conserved charges.")


def run(quick=False):
    masses=(.02,.025,.03,.04) if quick else (.02,.025,.03,.04,.05,.08,.1)
    results=[repaired_instantaneous_grid(x) for x in masses]
    return {
       "task":"FRW-indexed repair of spurious on-shell quantum collision invariants",
       "daughter_grid":"core 0..99 + coarse indices every 50 from 50 to selected upper cutoff, delta=.001",
       "instantaneous_mass_sweep":results,
       "finding":"A larger physical daughter UV domain can recover graph connectivity at larger M/T, but the minimum graph-valid cutoff does not give thermal-moment convergence.",
       "scope":"STATIC instantaneous graph search; no physical particle transport, no FRW energy closure, no CLASS."
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick",action="store_true")
    parser.add_argument("--output",default="outputs/dynamic_rank_uv_repair.json")
    opts=parser.parse_args()
    r=run(opts.quick)
    path=Path(opts.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(r,indent=2),encoding="utf-8")
    print(json.dumps(r,indent=2))
