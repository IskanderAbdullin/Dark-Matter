"""Explicit nonphysical MODULAR conservation law of a three-band kinetic grid.

Index set I = {0..15} U {2500..2519} U {5000..5019};
daughter energy e_i=(i+1/2)*delta and shared mother key k=i+j+1.
All daughter indices have residues r_i=i mod 2500 between 0 and 19,
so for every pair r_i+r_j<2500. Therefore
r_i+r_j = (i+j) mod 2500 = (k-1) mod 2500.

This produces an EXACT, nonphysical conserved pseudo-energy
Qmod = sum_i r_i N_phi_i + sum_j r_j N_chi2_j
       + sum_k ((k-1) mod 2500) N_chi1_k
in the full nonlinear semidiscrete collision equations, regardless of
occupations and coupling strength. It is independent of the three
physical energy/number conservation laws on this multi-band grid.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from collision_invariant_rank import three_block_indices
from adaptive_resonant_lattice import AdaptiveClosedDecay


def modular_invariant(model, period=2500):
    if period<2:
        raise ValueError("positive modular period >=2 required")
    index=model.index
    n=model.nq
    rb=index % period
    rf=rb.copy()
    ik,jk=np.nonzero(model.K)
    parent_keys=index[ik]+index[jk]+1
    unique=np.unique(parent_keys)
    expected=model.parent_E/model.delta
    # Expected parent-energy key must match integer resonance within fp.
    if np.max(abs(expected-np.rint(expected)))>1e-8:
        raise ArithmeticError("Mother energies are not exact integer multiples")
    mapped=np.rint(expected).astype(int)
    if not np.array_equal(mapped,unique):
        raise ArithmeticError("Mother resonance index mismatch")
    rp=(mapped-1)%period
    assert np.max(index%period)<period
    parent_idx=model.parent_index[ik,jk]
    residual=rb[ik]+rf[jk]-rp[parent_idx]
    max_err=float(np.max(abs(residual)))
    phys=np.zeros((2*n,3))
    phys[:n,0]=1.;phys[n:,1]=1.
    phys[:n,2]=model.e;phys[n:,2]=model.e
    q,_=np.linalg.qr(phys,mode="reduced")
    vector=np.r_[rb,rf].astype(float)
    independence=float(np.linalg.norm(vector-q@(q.T@vector)))
    return {
        "period":period,
        "daughter_residue_max":int(max(rb)),
        "every_active_reaction_modular_defect_max":max_err,
        "distance_from_span_of_two_charges_and_energy":independence,
        "boson_coefficients":rb,
        "chi2_coefficients":rf,
        "parent_coefficients":rp,
        "is_exact_extra_invariant":bool(max_err==0 and independence>1e-6),
    }


def check_full_nonlinear_rhs(model,weights):
    n=model.nq
    x=np.r_[0.1*np.exp(-model.e),
            0.2*np.exp(-model.e),
            1/(1+np.exp(model.parent_E))]
    dx=model.rhs(0.,x)
    return float(
      np.dot(weights["boson_coefficients"]*model.vol_phi,dx[:n])+
      np.dot(weights["chi2_coefficients"]*model.vol_chi,dx[n:2*n])+
      np.dot(weights["parent_coefficients"]*model.vol_parent,dx[2*n:]))


def run():
    model=AdaptiveClosedDecay(three_block_indices())
    result=modular_invariant(model)
    derivative=check_full_nonlinear_rhs(model,result)
    assert result["is_exact_extra_invariant"]
    assert abs(derivative)<1e-10
    return {
      "scope":"Exact algebraic witness for a grid-only, physically spurious collision invariant",
      "mesh":"{0..15} U {2500..2519} U {5000..5019}, daughter energies (i+.5)*delta, delta=.001",
      "reaction":"single mother state per unique i+j+1, chi1 <-> chi2+phi",
      "analytic_charge":"Q_mod=sum_i(i mod 2500)*N_phi_i + sum_j(j mod 2500)*N_chi2_j + sum_k((k-1) mod 2500)*N_chi1_k",
      "independent_from_two_number_charges_and_energy":result["is_exact_extra_invariant"],
      "max_reaction_modular_residual":result["every_active_reaction_modular_defect_max"],
      "nonlinear_ODE_modular_charge_derivative":derivative,
      "independence_norm":result["distance_from_span_of_two_charges_and_energy"],
      "physical_status":"NUMERICAL ARTIFACT; not a physical conservation law of the continuous Bose dark-matter model",
      "warning":"A mesh can satisfy the H theorem, correct detailed balance and physical invariants but still possess nonphysical conserved block-residue charge."
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/modular_collision_invariant.json")
    opts=p.parse_args()
    data=run()
    dst=Path(opts.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(data,indent=2),encoding="utf-8")
    print(json.dumps(data,indent=2))
