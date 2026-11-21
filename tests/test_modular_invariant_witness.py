"""A closed-form modular aliasing invariant must be caught by adaptive graph audits."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from modular_invariant_witness import modular_invariant,check_full_nonlinear_rhs
from adaptive_resonant_lattice import AdaptiveClosedDecay
from collision_invariant_rank import three_block_indices


def test_three_band_extra_charge_has_exact_modular_formula():
    m=AdaptiveClosedDecay(three_block_indices())
    w=modular_invariant(m,period=2500)
    assert w["is_exact_extra_invariant"]
    assert w["every_active_reaction_modular_defect_max"]==0
    assert w["daughter_residue_max"]==19
    assert w["distance_from_span_of_two_charges_and_energy"]>50
    assert abs(check_full_nonlinear_rhs(m,w))<1e-10


def test_modular_charge_conserved_for_nonthermal_random_population():
    m=AdaptiveClosedDecay(three_block_indices())
    w=modular_invariant(m,period=2500)
    rng=np.random.default_rng(1024)
    for _ in range(5):
        y=np.r_[rng.uniform(0.01,20,m.nq),
                rng.uniform(.01,.99,m.nq),
                rng.uniform(.01,.99,m.nparent)]
        d=m.rhs(0.,y)
        n=m.nq
        pseudo=float(np.dot(w["boson_coefficients"]*m.vol_phi,d[:n])+
                     np.dot(w["chi2_coefficients"]*m.vol_chi,d[n:2*n])+
                     np.dot(w["parent_coefficients"]*m.vol_parent,d[2*n:]))
        assert abs(pseudo)<1e-9


def test_modular_charge_is_not_universal_on_arbitrary_sparse_grids():
    from collision_invariant_rank import sparse_cold_hot_indices
    m=AdaptiveClosedDecay(sparse_cold_hot_indices())
    result=modular_invariant(m,period=2500)
    assert not result["is_exact_extra_invariant"]
    assert result["every_active_reaction_modular_defect_max"]>0
