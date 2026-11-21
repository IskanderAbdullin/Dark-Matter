"""Pinned BEST real RD-FRW time order: negative IR accuracy, positive work order."""
from __future__ import annotations
import json
import math
from pathlib import Path


def data():
    return json.loads((Path(__file__).resolve().parents[1]/
       "kinetics/results/best-rd-time-refinement-2026-11-21.json").read_text(encoding="utf-8"))


def test_external_best_4_8_16_32_rd_time_refinement_provenance():
    r=data()
    assert r["external_best_run"]==44919748455
    assert r["best_pin"]=="c70875be714e7a6bc6d38578c34bb0dc9f5aefa9"
    assert len(r["cases"])==5
    assert not r["claimed_time_convergence"]
    assert not r["full_86eV_cosmology_or_CLASS_ready"]
    assert math.isclose(r["final_a"],math.sqrt(1.4),rel_tol=1e-14)


def test_phi_number_ir_zero_mode_have_first_order_like_difference_ratios():
    comp=data()["half_step_comparisons"]
    assert [a["steps"] for a in comp]==["4->8","8->16","16->32"]
    for key in ("phi_n_relative_difference","qmin_f_relative_difference","N0_relative_difference"):
        mags=[abs(c[key]) for c in comp]
        assert mags[0]>mags[1]>mags[2]
        assert 1.9<mags[0]/mags[1]<2.1
        assert 1.9<mags[1]/mags[2]<2.1


def test_physical_pressure_work_trapezoid_is_second_order_not_exact():
    a=data()["half_step_comparisons"]
    assert all(3.9<x["pressure_work_residual_ratio"]<4.1 for x in a)
    assert .019<data()["deficits"]["soft_phi_qmin_16_to_32_rel"]<.021
    assert .0024<data()["deficits"]["phi_n_16_to_32_rel"]<.0027
    assert data()["interpreted_convergence_order"]["Richardson_is_not_physical_validation"]


def test_actual_best_energy_and_number_charge_residuals_are_separately_reported():
    for r in data()["cases"]:
        assert abs(r["charge1_relative_error"])<1e-8
        assert r["pressure_work_trapezoid"]>0
        assert math.isclose(
            r["physical_energy_change"]-r["pressure_work_trapezoid"],
            r["energy_ledger_residual"],rel_tol=1e-8,abs_tol=1e-13)
