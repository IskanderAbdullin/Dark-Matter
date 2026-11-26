"""Why the closed fixed-energy resonant 1<->2 stencil cannot be reused in FRW.

At reference x0=M1/T, each active reaction i+j<->k has
 Eparent_k(x0)=e_phi_i+e_chi_j. The freely streaming massive parent
corresponds to constant comoving momentum
 qparent_k = sqrt(Eparent_k(x0)^2-x0^2).

At later x, Eparent_k(x)=sqrt(qparent_k^2+x^2). For the SAME daughter
comoving momenta, Ephi_i+Echi_j is unchanged, so an old reaction has

 deltaE = sqrt(Eparent_k(x0)^2-x0^2+x^2)-Eparent_k(x0)
        = (x^2-x0^2)/(Eparent_k(x)+Eparent_k(x0)) != 0 for x!=x0.

This is an EXACT kinematic no-go for directly combining the
static-energy Schur closed solver with FRW redshift by merely
changing its rate prefactor: it would violate physical energy
conservation on every nonzero active channel.

The physical 1<->2 collision operator must use instantaneous parent
momentum/energy and reassign daughter energies, with a conservative
remapping/interpolation scheme that also preserves Pauli/Bose bounds,
two charges and an H theorem. This module only proves the obstruction.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import expit

from resonant_graph_repair import bridged_indices
from adaptive_resonant_lattice import AdaptiveClosedDecay


def resonance_drift(x0=.1,x=.2,max_index=5000):
    if x0<=0 or x<=0:
        raise ValueError("Positive mass/temperature ratios required")
    model=AdaptiveClosedDecay(
        bridged_indices(max_index=max_index),M1_T=x0)
    i,j=np.nonzero(model.K)
    k=model.parent_index[i,j]
    e0=model.parent_E
    if np.any(e0<x0*(1-1e-12)):
        raise ArithmeticError("Reference mother energy below mass")
    qparent=np.sqrt(np.maximum(e0*e0-x0*x0,0))
    e_new=np.hypot(qparent,x)
    daughter_e=model.e[i]+model.e[j]
    mismatch=e_new[k]-daughter_e
    algebra=(x*x-x0*x0)/(e_new[k]+e0[k])
    if np.max(abs(mismatch-algebra))>1e-12:
        raise ArithmeticError("Physical free-stream energy identity failed")
    rates=model.K[i,j]*expit(-e0[k])
    denominator=float(np.dot(rates,daughter_e))
    if denominator<=0:
        raise ArithmeticError("No nonzero active source rates")
    current_on_shell=4*model.e[i]*model.e[j]>=x*x
    return {
       "reference_x_M1_over_T":x0,
       "later_x_M1_over_T":x,
       "daughter_grid_max_index":max_index,
       "closed_graph_active_collisions":len(i),
       "shared_mother_states":model.nparent,
       "number_of_physical_collision_invariants_at_x0":
           model.graph["active_graph_collision_invariants"],
       "maximum_absolute_energy_defect_per_reference_channel":
           float(np.max(abs(mismatch))),
       "minimum_signed_energy_defect":float(np.min(mismatch)),
       "maximum_signed_energy_defect":float(np.max(mismatch)),
       "signed_weighted_energy_defect_relative_to_reference_source":
           float(np.dot(rates,mismatch)/denominator),
       "rate_weighted_fraction_newly_forbidden_by_instantaneous_kinematics":
           float(np.dot(rates,~current_on_shell)/rates.sum()),
       "unweighted_fraction_newly_forbidden_channels":
           float(np.mean(~current_on_shell)),
       "maximum_algebraic_energy_identity_residual":
           float(np.max(abs(mismatch-algebra))),
       "number_of_original_reactions_still_exactly_energy_resonant":
           int(np.count_nonzero(abs(mismatch)<1e-12)),
       "not_a_correct_FRW_collision_operator":True,
       "note":"Counterfactual reuse of the same reference energy triplets under collisionless mother redshift; this is NOT actual FRW production or a physical energy-loss prediction."
    }


def run():
    xs=(.1,.11,.2,.5,1.,2.)
    cases=[resonance_drift(x=x) for x in xs]
    assert cases[0]["maximum_absolute_energy_defect_per_reference_channel"]<1e-12
    for row in cases[1:]:
        assert row["minimum_signed_energy_defect"]>0,row
        assert row["number_of_original_reactions_still_exactly_energy_resonant"]==0,row
    return {
       "scope":"Exact kinematic obstruction to frozen resonant parent-energy grid under massive collisionless FRW redshift",
       "reference_mass_ratio_x0":.1,
       "cases":cases,
       "physical_conclusion":"A static set of daughter-energy and massive-parent-comoving momentum triplets cannot stay exactly resonant for more than one x=M1/T. Time-dependent on-shell remapping is mandatory.",
       "caveat":"Not a cosmological 1<->2 integration and not a primordial spectrum."
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/frw_resonance_drift_no_go.json")
    args=p.parse_args()
    data=run()
    path=Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2),encoding="utf-8")
    print(json.dumps(data,indent=2))
