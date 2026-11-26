"""Symmetric versus left-endpoint split-step FRW diagnostics."""
from __future__ import annotations
import sys
from pathlib import Path
import math
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from strang_frw_split import SymmetricFRWPilot,compare_lie_and_symmetric


def test_strang_substeps_conserve_both_charges_and_discrete_collision_energy():
    result=SymmetricFRWPilot(steps=4).run()
    assert result["global_two_charge_drift"]<1e-8
    assert result["final"]["n_phi"]>0
    assert result["CLASS_ready"] is False
    assert result["physical_total_energy_conservation_claimed"] is False
    assert len(result["steps"])==4
    for row in result["steps"]:
        assert row["maximum_charge_error"]<1e-8
        assert abs(row["left_collision_discrete_energy_error"])<1e-8
        assert abs(row["right_collision_discrete_energy_error"])<1e-8
        assert math.isfinite(row["left_collision_physical_shell_error"])
        assert math.isfinite(row["right_collision_physical_shell_error"])
        assert math.isfinite(row["numerical_remap_energy_error"])


def test_symmetric_scheme_allows_rebuilding_parent_basis_with_fixed_q_domain():
    result=SymmetricFRWPilot(x0=.02,xend=.022,steps=4,kappa=10.,
                              allow_mother_mode_changes=True).run()
    assert result["parent_modes_final"]<result["parent_modes_initial"]
    assert result["global_two_charge_drift"]<1e-8
    assert result["final"]["n_chi2"]>0
    assert result["final"]["n_phi"]>0
    assert result["time_refinement_not_assumed"]


def test_lie_and_strang_refinement_are_reported_without_unsupported_order_claim():
    result=compare_lie_and_symmetric()
    assert result["no_physical_time_convergence_claim"]
    assert len(result["left_endpoint_Lie"]["runs"])==2
    assert len(result["symmetric"]["runs"])==2
    for method in ("left_endpoint_Lie","symmetric"):
        d=result[method]["relative_step_refinement_difference"]
        assert math.isfinite(d)
        assert 0<=d<2
        assert all(row["charge_drift"]<1e-8 for row in result[method]["runs"])
    assert isinstance(result["symmetric_has_smaller_observed_step_error"],bool)
