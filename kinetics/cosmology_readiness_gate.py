"""Fail-closed readiness gate for handing a kinetic spectrum to CLASS.

Numerical preservation of energy, charge, detailed balance and the H
theorem does NOT imply physical correctness. Especially, a sampled
reaction network must have no unintended collision invariants.
All independently required physics/convergence evidence must be supplied
explicitly; absent checks are NOT silently treated as passing.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from collision_invariant_rank import (
    audit_integer_grid, sparse_cold_hot_indices, three_block_indices
)


REQUIRED_EXTERNAL_CHECKS=(
    "full_spectral_convergence",
    "infrared_regulator_and_zero_mode_control",
    "closed_FRW_energy_and_entropy",
    "realistic_mother_depletion_and_rethermalization",
    "physical_masses_and_width",
    "correct_relic_density_without_rescaling",
)


def evaluate_candidate(graph, parent_inventory_relative_bias,
                       evidence=None, parent_bias_tolerance=.01):
    if not isinstance(graph,dict):
        raise TypeError("Expected collision-graph audit dictionary")
    evidence=evidence or {}
    if parent_bias_tolerance<=0:
        raise ValueError("positive tolerance required")
    flags={
       "collision_graph_has_exactly_three_invariants":
           graph.get("additional_invariants_beyond_number_number_energy")==0,
       "no_inactive_mother_modes":
           graph.get("inactive_declared_parent_modes")==0,
       "no_isolated_daughter_modes":
           graph.get("isolated_boson_modes")==0 and graph.get("isolated_chi2_modes")==0,
       "parent_thermal_inventory_within_tolerance":
           parent_inventory_relative_bias is not None and
           abs(parent_inventory_relative_bias)<parent_bias_tolerance,
    }
    flags.update({key:evidence.get(key) is True for key in REQUIRED_EXTERNAL_CHECKS})
    reasons=[key for key,passed in flags.items() if not passed]
    return {
        "approved_for_CLASS":len(reasons)==0,
        "passed":flags,
        "failed_or_unverified":reasons,
        "numerical_parent_inventory_bias":parent_inventory_relative_bias,
        "parent_bias_tolerance":parent_bias_tolerance,
        "warning":"A pass only confirms supplied gates, not observational viability. Never invent evidence for unverified checks."
    }


def require_ready_for_CLASS(report):
    """Raise to prohibit accidentally exporting an unvalidated f(q)."""
    if not report.get("approved_for_CLASS",False):
        raise ValueError("CLASS initial PSD rejected: "+", ".join(
            report.get("failed_or_unverified",["no positive evidence supplied"])))
    return True


def run():
    experiments={
      "sparse_cold_hot_27":(audit_integer_grid(sparse_cold_hot_indices()),.0023424248092460687),
      "three_blocks_56":(audit_integer_grid(three_block_indices()),-.20621124318297168)
    }
    out={}
    for name,(graph,bias) in experiments.items():
        out[name]=evaluate_candidate(graph,bias)
        assert not out[name]["approved_for_CLASS"]
    return {
      "scope":"Explicit fail-closed regression gate for primordial PSD release",
      "evidence_policy":"Missing evidence is rejected; script cannot certify the actual physics",
      "required_external_checks":list(REQUIRED_EXTERNAL_CHECKS),
      "candidate_results":out
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/kinetic_class_readiness.json")
    opts=p.parse_args()
    result=run()
    dst=Path(opts.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
