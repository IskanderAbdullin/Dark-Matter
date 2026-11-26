"""Exact collision-invariant audit: physical vs grid-induced conserved charges."""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from collision_invariant_rank import (
    audit_integer_grid, three_block_indices, sparse_cold_hot_indices
)


@pytest.mark.parametrize("n,delta,unreachable",[(30,.08,1),(50,.08,1),(30,.02,4)])
def test_uniform_graph_has_only_three_invariants_after_unused_parent_removal(n,delta,unreachable):
    result=audit_integer_grid(range(n),delta=delta,
                              declared_parent_keys=range(1,2*n))
    assert result["active_graph_collision_invariants"]==3
    assert result["inactive_declared_parent_modes"]==unreachable
    assert result["full_grid_collision_invariants_including_inactive_parents"]==3+unreachable
    assert result["analytical_three_invariant_residual"]<1e-12


def test_cold_sparse_grid_has_eleven_extra_collision_invariants():
    x=audit_integer_grid(sparse_cold_hot_indices())
    assert x["daughter_grid_n"]==27
    assert x["active_graph_collision_invariants"]==14
    assert x["additional_invariants_beyond_number_number_energy"]==11
    assert x["isolated_boson_modes"]==0
    assert x["isolated_chi2_modes"]==0
    assert not x["accept_for_unique_chemical_equilibrium"]
    witness=x["extra_invariant_witness"]
    assert witness is not None
    assert witness["extra_invariant_max_reaction_constraint_residual"]<1e-8
    assert witness["orthogonality_to_three_physical_invariants"]<1e-8


def test_three_block_grid_has_one_extra_invariant_even_with_all_daughters_active():
    x=audit_integer_grid(three_block_indices())
    assert x["daughter_grid_n"]==56
    assert x["active_graph_collision_invariants"]==4
    assert x["additional_invariants_beyond_number_number_energy"]==1
    assert x["inactive_declared_parent_modes"]==0
    assert x["isolated_boson_modes"]==0
    assert x["isolated_chi2_modes"]==0
    assert x["analytical_three_invariant_residual"]<1e-12


def test_invariants_are_robust_to_scale_of_energy_units_for_same_graph():
    a=audit_integer_grid(three_block_indices(),delta=.001,mass=.1)
    b=audit_integer_grid(three_block_indices(),delta=.01,mass=1.)
    assert a["active_pair_reactions"]==b["active_pair_reactions"]
    assert a["active_graph_collision_invariants"]==b["active_graph_collision_invariants"]


def test_production_gate_rejects_inactive_parent_modes():
    x=audit_integer_grid(range(30),delta=.08,
                         declared_parent_keys=range(1,60))
    assert x["additional_invariants_beyond_number_number_energy"]==0
    assert x["inactive_declared_parent_keys"]==[1]
    assert not x["accept_for_unique_chemical_equilibrium"]
