"""Actual independent 64/128 step BEST FRW Bose data and analytic redshift audit."""
from __future__ import annotations
import json
import math
from pathlib import Path


def read():
    path=(Path(__file__).resolve().parents[1]/
          "kinetics/results/best-rd-time-64-128-real-2026-11-21.json")
    return json.loads(path.read_text(encoding="utf-8"))


def test_real_independent_best_evidence_is_pinned_and_not_class_ready():
    d=read()
    assert d["actual_BEST_run_id"]==44923390224
    assert d["external_git_sha"]=="c70875be714e7a6bc6d38578c34bb0dc9f5aefa9"
    assert d["actual_mpi_ranks"]==2
    assert d["raw_zip_artifact_id"]==14714124787
    assert d["claims"]["true_external_BEST_time_evolution"]
    assert not d["claims"]["fully_time_converged"]
    assert not d["claims"]["macroscopic_BEC_established"]
    assert not d["claims"]["primordial_86eV_dark_matter_CLASS_ready"]


def test_64_to_128_time_refinement_remains_first_order_in_soft_bose_f():
    d=read()
    c64,c128=d["cases"]
    assert c64["steps"]==64 and c128["steps"]==128
    assert math.isclose(c64["dtGamma"],2*c128["dtGamma"])
    assert math.isclose(d["final_a"],math.sqrt(1.4),rel_tol=1e-14)
    delta=d["comparisons"]
    assert math.isclose(delta["phi_soft_relative_64_to_128"],
                        c64["phi_f_at_qmin"]/c128["phi_f_at_qmin"]-1,
                        rel_tol=1e-12)
    assert .0049<delta["phi_soft_relative_64_to_128"]<.0052
    assert -.0007<delta["phi_number_relative_64_to_128"]<-.0006


def test_independent_comoving_parent_analytic_limit_predicts_BEST_extrapolation():
    c=read()["comparisons"]
    analytic=c["analytic_f_qmin"]
    extrap=c["first_order_Richardson_64_128_f"]
    assert math.isclose(extrap,1.019468120189062,rel_tol=1e-12)
    assert abs(extrap/analytic-1)<.0003
    assert .004<c["best_128_f_to_analytic_relative"]<.006


def test_global_conservation_does_not_invent_ir_accuracy():
    r=read()
    for row in r["cases"]:
        assert abs(row["relative_chi1_phi_charge_error"])<1e-9
        assert abs(row["relative_chi1_chi2_charge_error"])<1e-9
        assert math.isclose(
            row["energy_change"]-row["pressure_work_trapezoid"],
            row["energy_work_residual"],rel_tol=1e-7,abs_tol=1e-13)
    assert r["comparisons"]["energy_work_residual_64_to_128_ratio"]>3
    assert not r["claims"]["fully_time_converged"]
