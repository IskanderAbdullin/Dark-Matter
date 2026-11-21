"""Reject fake/missing independent references and test source-stoichiometry math."""
from __future__ import annotations
import copy
import math
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from best_external_vacuum_reference import BEST_PIN,benchmark_protocol
from best_source_comparison import inspect_record,evaluate_records


def synthetic_fixture(grid_n=100,neval=20000,phi_fractional_error=0.):
    """TEST-ONLY fabricated input — absolutely no BEST execution."""
    p=benchmark_protocol(grid_n=grid_n,vegas_neval=neval)
    target=p["independently_derived_vacuum_source_n_over_T3_per_time"]
    rate=target*(1+phi_fractional_error)
    return {
        "external_solver_executed":True,
        "BEST_source_commit_verified":BEST_PIN,
        "BEST_git_revision_matched_pinned_protocol":True,
        "model_and_protocol":p,
        "instantaneous_integrated_number_rates":{
            "chi1":-rate,"phi":rate,"chi2":rate}
    }


def test_no_real_external_files_means_unverified():
    result=evaluate_records([])
    assert result["source_agreement_verified"] is False
    assert result["records_checked"]==0
    assert "UNVERIFIED" in result["status"]
    assert result["freezeout_or_CLASS_validated"] is False


def test_single_perfect_fake_math_fixture_is_not_independent_validation():
    case=synthetic_fixture()
    result=evaluate_records([case])
    assert result["checks"][0]["single_record_passes_preliminary_criteria"]
    assert result["source_agreement_verified"] is False
    assert "UNVERIFIED" in result["status"]


def test_four_synthetic_records_test_math_gate_only_not_true_provenance():
    inputs=[synthetic_fixture(nq,neval)
            for nq in (80,120) for neval in (20000,50000)]
    report=evaluate_records(inputs)
    assert report["source_agreement_verified"]
    assert report["freezeout_or_CLASS_validated"] is False
    assert report["multiple_independent_momentum_grids"]
    assert report["multiple_monte_carlo_sample_counts"]
    assert len(report["checks"])==4
    # These are deliberately fabricated unit fixtures, not BEST runs.


def test_incorrect_signs_and_missing_chi2_rate_are_not_accepted():
    a=synthetic_fixture()
    a["instantaneous_integrated_number_rates"]["chi1"]*= -1
    r=inspect_record(a)
    assert r["correct_physical_signs"] is False
    assert not r["single_record_passes_preliminary_criteria"]
    b=synthetic_fixture()
    del b["instantaneous_integrated_number_rates"]["chi2"]
    with pytest.raises(ValueError,match="Missing"):
        inspect_record(b)


def test_bad_source_rate_stoichiometry_and_unpinned_commit_are_rejected():
    a=synthetic_fixture(phi_fractional_error=.3)
    r=inspect_record(a)
    assert not r["initial_source_within_rate_tolerance"]
    b=synthetic_fixture()
    b["instantaneous_integrated_number_rates"]["chi2"]*=.7
    assert not inspect_record(b)["all_three_numbers_within_charge_tolerance"]
    c=synthetic_fixture()
    c["BEST_source_commit_verified"]="0"*40
    with pytest.raises(ValueError,match="pinned"):
        inspect_record(c)
    d=synthetic_fixture()
    d["external_solver_executed"]=False
    with pytest.raises(ValueError,match="Unexecuted"):
        inspect_record(d)


def test_spin_summed_matrix_element_mismatch_fail_closed():
    p=synthetic_fixture()
    p["model_and_protocol"]["best_spin_summed_matrix_element_squared"]*=.5
    with pytest.raises(ValueError,match="matrix-element"):
        inspect_record(p)


def test_nonfinite_and_negative_tolerance_fail_closed():
    a=synthetic_fixture()
    a["instantaneous_integrated_number_rates"]["phi"]=float("nan")
    with pytest.raises(ValueError,match="Non-finite"):
        inspect_record(a)
    with pytest.raises(ValueError):
        inspect_record(synthetic_fixture(),rate_tolerance=0.)
