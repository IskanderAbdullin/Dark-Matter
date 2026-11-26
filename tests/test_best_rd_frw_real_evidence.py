"""Actual prescribed-radiation BEST source evidence and negative dt-convergence."""
from __future__ import annotations
import math
import json
from pathlib import Path


def data():
    return json.loads((Path(__file__).resolve().parents[1]/
        "kinetics/results/best-rd-frw-real-2026-11-21.json").read_text(encoding="utf-8"))


def test_external_best_real_expanding_frw_run_and_physical_scope():
    d=data()
    assert d["actual_independent_best_executed"]
    assert d["run_id"]==44919387818
    assert d["best_sha"]=="c70875be714e7a6bc6d38578c34bb0dc9f5aefa9"
    assert math.isclose(d["final_a"],math.sqrt(1.4),rel_tol=1e-13)
    assert d["time_step_convergence_claimed"] is False
    assert d["physical_86eV_relic_or_CLASS_validated"] is False
    assert not d["SM_plasma_or_gstarS_change_included"]


def test_both_charges_and_physical_energy_work_are_reported_not_hidden():
    d=data()
    for c in d["cases"]:
        assert abs(c["relative_charge1_error"])<1e-8
        assert abs(c["relative_charge2_error"])<1e-8
        assert c["physical_RD_pressure_work_approx"]>0
        assert math.isclose(
            c["comoving_energy_change"]-c["physical_RD_pressure_work_approx"],
            c["physical_energy_residual"],rel_tol=1e-8,abs_tol=1e-13)
        assert math.isclose(c["relative_energy_residual"],
                            c["physical_energy_residual"]/d["comoving_energy_initial"],
                            rel_tol=1e-7)


def test_dt_refinement_improves_pressure_work_but_weakly_converges_cold_source():
    d=data()
    errors=d["errors"]
    assert .009<abs(errors["exprb_4_to_8_phi_relative_difference"])<.012
    assert .07<errors["exprb_4_to_8_qmin_f_relative_difference"]<.09
    assert .07<errors["exprb_4_to_8_zero_mode_N0_relative_difference"]<.09
    assert 3.8<errors["physical_work_quadrature_residual_factor_improvement"]<4.2
    assert not d["physical_FRW_energy_balance_validated_to_asymptotic_limit"]


def test_closed_and_expanding_bose_source_differs_from_fixed_a_without_rethermalizing():
    d=data()
    c=d["cases"][1]
    assert c["final_n_phi"]>7e-7
    assert c["final_phi_f_qmin"]>1
    assert c["final_phi_n0"]>1e-14
    assert not d["zero_mode_macroscopic_BEC_confirmed"]
