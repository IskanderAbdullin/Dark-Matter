"""Two independent angular parameterizations of identical decay kinematics.

This scans the entropy-balanced SURROGATE, never the original physical ODE.
No FRW, finite-mass daughters, mother depletion, condensate or CLASS.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from lattice_detailed_balance import EntropyBalancedLattice


def one(nq,np1,nmu,rule):
    model=EntropyBalancedLattice(nq,np1,nmu,angular_rule=rule)
    return {
      "grid":[nq,np1,nmu],"rule":rule,
      "min_generated_daughter_q":model.min_quadrature_daughter_q,
      "max_angular_normalization_error":model.angular_normalization_error,
      "original_parent_vacuum_rate":model.true_vacuum_rate,
      "surrogate_vacuum_rate":model.lattice_vacuum_rate,
      "physical_parent_energy_projection_error":model.absolute_parent_energy_error,
      "structural_checks":model.checks(),
      "evolution":model.evolve()}


def run(quick=False):
    cases=([(100,140,24),(100,140,48)] if quick else
          [(100,240,24),(100,240,48),(100,240,96),
           (150,240,24),(150,240,48),(150,240,96),
           (220,240,24),(220,240,48),(220,240,72),
           (220,240,96),(220,240,120),(220,240,160),
           (220,240,240),(300,240,48),(300,240,96)])
    rows=[one(*case,rule) for case in cases for rule in ("mu","logq")]
    for r in rows:
        c=r["structural_checks"]
        assert c["thermal_equilibrium_max_bin_rate"]<1e-11
        assert abs(c["daughter_number_rate_mismatch"])<1e-10
    return {
      "model":"fixed-bath 1<->2 lattice surrogate, comparative angular quadrature",
      "reference_commit":"ea6ec461997600d5f00bf244ea808607d40ef8dd",
      "scope_warning":"Mother thermal occupation uses effective pair energy; this is NOT a physically exact collision integral.",
      "main_test":"At fixed q bin resolution and parent quadrature, log-q integration resolves backward collinear emission without exponentially increasing angular nodes.",
      "results":rows
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick",action="store_true")
    parser.add_argument("--output",default="outputs/logq_angular_convergence.json")
    opts=parser.parse_args()
    result=run(opts.quick)
    p=Path(opts.output);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
