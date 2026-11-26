"""Narrow early-FRW split: rebuild instantaneous graph and audit all budgets."""
from __future__ import annotations
import sys
from pathlib import Path
import math
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from narrow_closed_frw_split import NarrowClosedFRWPilot


def test_instantaneous_on_shell_network_changes_with_temperature():
    p=NarrowClosedFRWPilot(steps=4)
    assert p.models[0].graph["active_graph_collision_invariants"]==3
    assert p.models[-1].graph["active_graph_collision_invariants"]==3
    assert p.models[-1].active_pairs<p.models[0].active_pairs
    assert all(m.max_resonance_error<1e-12 for m in p.models)
    assert all(np.array_equal(m.parent_E,p.models[0].parent_E) for m in p.models)


def test_four_steps_evolve_3_species_and_keep_comoving_charges():
    p=NarrowClosedFRWPilot(steps=4)
    r=p.run()
    assert r["not_a_self_consistent_cosmological_model"]
    assert r["final"]["n_phi"]>0
    assert r["final"]["n_chi2"]>0
    assert r["shared_parent_modes"]>300
    assert r["max_global_charge_drift"]<1e-8
    assert len(r["steps"])==4
    for row in r["steps"]:
        assert row["maximum_charge_drift"]<1e-8
        assert abs(row["collision_discrete_energy_error"])<1e-8
        assert np.isfinite(row["collision_physical_shell_energy_error"])
        assert np.isfinite(row["collisionless_FRW_energy_change"])
        assert np.isfinite(row["finite_volume_remap_energy_error"])
        assert row["min_fermion_occupation"]>=-1e-8
        assert row["max_fermion_occupation"]<=1+1e-8
    assert r["maximum_energy_remap_defect"]>=0


def test_step_refinement_keeps_total_nphi_approximately_consistent():
    coarse=NarrowClosedFRWPilot(steps=2).run()
    fine=NarrowClosedFRWPilot(steps=4).run()
    n1=coarse["final"]["n_phi"];n2=fine["final"]["n_phi"]
    assert n1>0 and n2>0
    # The preliminary Lie split is NOT converged: the 2-to-4 step
    # change is ~25.7% for these parameters. Record its failure instead
    # of silently asserting a physically unjustified 5% accuracy.
    drift=abs(n1/n2-1)
    assert .1<drift<.4
    assert fine["time_step_convergence_not_established"]
    assert fine["max_global_charge_drift"]<1e-8


def test_rejects_invalid_or_threshold_crossing_parent_energy_grids():
    with pytest.raises(ValueError):
        NarrowClosedFRWPilot(x0=.02,xend=.02)
    with pytest.raises(ValueError):
        NarrowClosedFRWPilot(x0=.02,xend=.025)
    with pytest.raises(ValueError):
        NarrowClosedFRWPilot(kappa=0)
