"""Real 2-stage BEST exprb2 FRW microphysics; numerical order not relic claim."""
from __future__ import annotations
import json
import math
from pathlib import Path


def data():
    f=Path(__file__).resolve().parents[1]/"kinetics/results/best-frw-exprb2-real-2026-11-21.json"
    return json.loads(f.read_text(encoding="utf-8"))


def test_real_external_best_second_order_source_has_verified_provenance():
    d=data()
    assert d["run_id"]==44924314432
    assert d["actual_best_git_sha"]=="c70875be714e7a6bc6d38578c34bb0dc9f5aefa9"
    assert d["independent_mpi_ranks"]==2
    assert len(d["cases"])==4
    assert [x["steps"] for x in d["cases"]]==[4,8,16,32]
    assert d["integrator"].startswith("BEST v1.3.5 native exprb2")
    assert not d["frozen_86eV_cosmological_DM_and_CLASS_ready"]
    assert not d["macroscopic_BEC_claimed"]


def test_soft_bose_first_vs_second_order_error_character_is_observed():
    d=data()
    differences=[abs(x["f_qmin_rel"]) for x in d["comparisons"]]
    assert differences[0]>differences[1]>differences[2]
    # With Monte Carlo noise and finite q grid, the second-order pattern
    # is approximate (not an enforced perfect ratio of four).
    assert 3.5<differences[0]/differences[1]<6.0
    assert 3.5<differences[1]/differences[2]<6.0
    assert all(3.8<x["pressure_work_residual_ratio"]<4.2
               for x in d["comparisons"])


def test_32_exprb2_steps_match_independent_analytic_redshift_gain():
    d=data()
    c32=d["cases"][-1]
    ref=d["independent_FRW_collisionless_mother_analytic_f_qmin"]
    assert math.isclose(ref,1.0196839876614525,rel_tol=1e-11)
    assert math.isclose(c32["final_phi_qmin_f"],1.019622172996092,rel_tol=1e-11)
    assert abs(c32["final_phi_qmin_f"]/ref-1)<1e-4
    assert abs(d["BEST_exprb2_32_N0_relative_to_independent_regular_cell"])<.02
    assert abs(c32["relative_charge1_residual"])<1e-8


def test_second_order_32_steps_beats_first_order_128_on_soft_source():
    d=data()
    reference=d["independent_FRW_collisionless_mother_analytic_f_qmin"]
    newer=d["cases"][-1]["final_phi_qmin_f"]
    old=d["previous_128_step_FIRST_ORDER_exprb_phi_qmin_f"]
    assert abs(newer/reference-1)<abs(old/reference-1)/40
    assert not d["qmin_kinetic_FULL_CONVERGENCE_established"]
    assert not d["SM_thermalization_and_heat_included"]
