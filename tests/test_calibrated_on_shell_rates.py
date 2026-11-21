"""Positive thermal cut-cell calibration without unphysical new channels."""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.special import expit

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from calibrated_on_shell_rates import ThermalCutCellCalibratedDecay,run
from resonant_graph_repair import bridged_indices
from sparse_kinetic_jacobian import SparseClosedDecay
from schur_closed_bdf import solve_closed_schur


@pytest.mark.parametrize("maximum",[5000,8500])
def test_calibrated_source_equals_high_order_FD_integral_over_retained_channels(maximum):
    indices=bridged_indices(max_index=maximum)
    cal=ThermalCutCellCalibratedDecay(indices)
    old=SparseClosedDecay(indices)
    assert cal.graph["active_graph_collision_invariants"]==3
    assert cal.active_pairs==old.active_pairs
    assert np.array_equal(cal.i,old.i) and np.array_equal(cal.j,old.j)
    assert np.array_equal(cal.k,old.k)
    assert np.all(cal.w>0)
    assert cal.reference_thermal_source_on_allowed_centroids>0
    assert abs(cal.initial_source_calibrated/
               cal.reference_thermal_source_on_allowed_centroids-1)<1e-12
    assert cal.excluded_centroid_cutcell_source>=0


def test_calibrated_ODE_still_has_conserved_charges_energy_and_h_theorem():
    model=ThermalCutCellCalibratedDecay(bridged_indices(max_index=5000))
    checks=model.thermodynamic_checks()
    assert checks["equilibrium_residual"]<1e-11
    assert checks["number_charge_derivative_max_abs"]<1e-10
    assert checks["energy_derivative_abs"]<1e-10
    assert checks["negative_entropy_rate"]<0
    assert checks["entropy_production"]>0
    assert abs(checks["entropy_balance_residual"])<1e-8


def test_sparse_jacobian_matches_finite_difference_after_rate_calibration():
    model=ThermalCutCellCalibratedDecay(bridged_indices(max_index=5000))
    state=np.r_[.02*np.exp(-model.e),.1*np.exp(-model.e),
                expit(-model.parent_E)]
    v=np.random.default_rng(2026).normal(size=model.size)
    assert model.jacobian_directional_error(state,v)<1e-6
    assert np.max(abs(model.rhs(0,state)-model.rhs_sparse(0,state)))<1e-7


def test_closed_dynamic_calibration_changes_kinetics_but_not_invariants():
    model=ThermalCutCellCalibratedDecay(bridged_indices(max_index=5000))
    out=solve_closed_schur(model,tend=.0002)
    assert out["graph_collision_invariants"]==3
    assert out["n_phi_over_T3"]>0
    assert out["n_mother_over_T3"]<model.mother_thermal_number_discretized
    assert max(abs(v) for v in out["number_and_energy_residuals"].values())<1e-8
    assert out["min_occupation"]>=-1e-7


def test_quick_experiment_does_not_claim_complete_physical_model():
    out=run(quick=True)
    assert not out["not_a_frozen_dark_matter_PSD"] is False
    assert out["graph_changed"] is False
    assert len(out["cases"])==1
    assert "thermal source" in out["caveat"]
