"""Independent actual 3-species massive BEST time evolution with corrected seed."""
from __future__ import annotations
import json
import math
from pathlib import Path


def data():
    return json.loads((Path(__file__).resolve().parents[1]/
        "kinetics/results/best-closed-massive-dynamic-real-2026-11-21.json").read_text(
        encoding="utf-8"))


def test_actual_pinned_best_mpi_dynamic_run_and_frw_class_disclaimers():
    d=data()
    assert d["run_id"]==44918653040
    assert d["real_independent_best_v1_3_5_sha"]=="c70875be714e7a6bc6d38578c34bb0dc9f5aefa9"
    assert d["actual_Mpi_ranks"]==2
    assert d["scoped_claims"]["real_external_kinetic_time_evolution"]
    assert d["scoped_claims"]["dynamic_massive_bose_zero_mode_N0_reported"]
    assert d["scoped_claims"]["FRW_expansion_or_SM_heat_exchange_included"] is False
    assert d["scoped_claims"]["true_macroscopic_BEC_established"] is False
    assert d["scoped_claims"]["relic_86eV_or_CLASS_established"] is False
    assert "exp(-q)" in d["initial_occupation_seed"]


def test_negligible_physical_initial_daughter_population_after_seed_fix():
    d=data()
    initial=d["initial_numbers"]
    assert initial["phi"]<2e-13
    assert initial["chi2"]<5e-13
    assert initial["chi1"]>.36
    assert math.isclose(d["initial_total_energy"],1.1506398648173843,rel_tol=1e-15)


def test_exponential_bose_gain_agrees_better_than_half_percent_with_full_closed_best():
    d=data()
    r=d["analytic_open_reservoir_early"]
    assert math.isclose(r["source_A_qmin"],27046.973477458418,rel_tol=1e-12)
    expected_f=r["predicted_exponential_f_qmin"]
    expected_N0=r["predicted_below_qmin_number"]
    for case in d["experiments"]:
        assert abs(case["final_qmin_f"]/expected_f-1)<.002
        assert abs(case["final_n0"]/expected_N0-1)<.02
        assert case["final_phi_n"]>0
        assert case["final_qmin_f"]>1
        assert case["final_n0"]<1e-12


def test_both_closed_charges_and_physical_energy_have_small_relative_drift():
    for r in data()["experiments"]:
        assert abs(r["charge_phi_rel"])<1e-8
        assert abs(r["charge_chi2_rel"])<1e-8
        assert abs(r["energy_rel"])<1e-9


def test_time_step_refinement_difference_and_separate_monte_carlo_caution():
    d=data()
    a,b=d["experiments"][:2]
    assert a["steps"]==4 and b["steps"]==8
    assert abs(a["final_phi_n"]/b["final_phi_n"]-1)<1e-4
    assert abs(a["final_qmin_f"]/b["final_qmin_f"]-1)<.005
    assert d["scoped_claims"]["time_qmin_qmax_seed_convergence_established"] is False
    assert d["diagnostics"]["actual_closed_phi_n_more_than_vacuum_linear"]>0
