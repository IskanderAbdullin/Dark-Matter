"""No unvalidated Bose kinetic spectrum can silently become CLASS input."""
from __future__ import annotations
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from collision_invariant_rank import audit_integer_grid,three_block_indices
from cosmology_readiness_gate import evaluate_candidate,require_ready_for_CLASS


def test_sparse_mesh_is_rejected_even_with_perfect_mother_normalization():
    graph=audit_integer_grid(three_block_indices())
    evidence={key:True for key in (
        "full_spectral_convergence",
        "infrared_regulator_and_zero_mode_control",
        "closed_FRW_energy_and_entropy",
        "realistic_mother_depletion_and_rethermalization",
        "physical_masses_and_width",
        "correct_relic_density_without_rescaling"
    )}
    report=evaluate_candidate(graph,0.0,evidence=evidence)
    assert not report["approved_for_CLASS"]
    assert "collision_graph_has_exactly_three_invariants" in report["failed_or_unverified"]
    with pytest.raises(ValueError,match="CLASS initial PSD rejected"):
        require_ready_for_CLASS(report)


def test_uniform_active_graph_still_rejected_without_physical_evidence():
    graph=audit_integer_grid(range(30),delta=.08)
    report=evaluate_candidate(graph,0.001)
    assert report["passed"]["collision_graph_has_exactly_three_invariants"]
    assert report["passed"]["parent_thermal_inventory_within_tolerance"]
    assert not report["approved_for_CLASS"]
    assert "full_spectral_convergence" in report["failed_or_unverified"]
    with pytest.raises(ValueError):
        require_ready_for_CLASS(report)


def test_missing_studies_never_become_implicit_success():
    graph=audit_integer_grid(range(30),delta=.08)
    report=evaluate_candidate(graph,None,evidence={"full_spectral_convergence":True})
    assert "parent_thermal_inventory_within_tolerance" in report["failed_or_unverified"]
    assert not report["approved_for_CLASS"]


def test_explicit_all_pass_is_possible_but_not_observational_validation():
    graph=audit_integer_grid(range(30),delta=.08)
    from cosmology_readiness_gate import REQUIRED_EXTERNAL_CHECKS
    report=evaluate_candidate(graph,0.0,{key:True for key in REQUIRED_EXTERNAL_CHECKS})
    assert report["approved_for_CLASS"]
    assert require_ready_for_CLASS(report) is True


def test_nonpositive_parent_tolerance_rejected():
    with pytest.raises(ValueError):
        evaluate_candidate(audit_integer_grid(range(30),delta=.08),0.0,parent_bias_tolerance=0)



def test_physically_unimplemented_plasma_and_late_time_checks_block_CLASS_even_if_old_checks_pass():
    graph=audit_integer_grid(range(30),delta=.08)
    from cosmology_readiness_gate import REQUIRED_EXTERNAL_CHECKS
    historical_checks=(
        "full_spectral_convergence",
        "infrared_regulator_and_zero_mode_control",
        "closed_FRW_energy_and_entropy",
        "realistic_mother_depletion_and_rethermalization",
        "physical_masses_and_width",
        "correct_relic_density_without_rescaling"
    )
    evidence={key:True for key in historical_checks}
    assert set(historical_checks)<set(REQUIRED_EXTERNAL_CHECKS)
    report=evaluate_candidate(graph,0.0,evidence)
    assert not report["approved_for_CLASS"]
    assert "physical_SM_heat_reservoir_and_entropy_gstarS" in report["failed_or_unverified"]
    assert "microphysical_bath_rates_and_stimulated_2to2_audit" in report["failed_or_unverified"]
    assert "late_time_hot_tail_and_BEC_state_convergence" in report["failed_or_unverified"]
    with pytest.raises(ValueError):
        require_ready_for_CLASS(report)


def test_partial_sm_bath_pilot_and_no_lagrangian_is_still_rejected():
    graph=audit_integer_grid(range(30),delta=.08)
    from cosmology_readiness_gate import REQUIRED_EXTERNAL_CHECKS
    evidence={key:True for key in REQUIRED_EXTERNAL_CHECKS}
    evidence["physical_SM_heat_reservoir_and_entropy_gstarS"]=False
    evidence["microphysical_bath_rates_and_stimulated_2to2_audit"]=False
    report=evaluate_candidate(graph,0,evidence)
    assert not report["approved_for_CLASS"]
