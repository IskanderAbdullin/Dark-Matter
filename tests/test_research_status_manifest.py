"""Research evidence manifest must retain blockers and negative numerical findings."""
from __future__ import annotations
import json
from pathlib import Path


def manifest():
    p=Path(__file__).resolve().parents[1]/"kinetics"/"results"/"research-status-2026-11-21.json"
    return json.loads(p.read_text(encoding="utf-8"))


def test_no_observational_claims_without_real_external_cross_validation():
    m=manifest()
    assert m["CLASS_initial_spectrum_approved"] is False
    assert m["diagnostics"]["best_external_executed"] is False
    assert m["diagnostics"]["best_independent_agreement_confirmed"] is False
    assert "NOT EXECUTED" in m["external_BEST_runtime_status"]
    assert any("CLASS" in s for s in m["disallowed_claims"])


def test_negative_frw_time_step_result_and_recovery_commit_are_preserved():
    m=manifest()
    assert .25<m["diagnostics"]["Lie_relative_drift"]<.27
    steps=m["sequential_commits_this_continuation"]
    assert len(steps)==14
    assert len({r["sha"] for r in steps})==14
    assert any(r["sha"].startswith("24cba826") for r in steps)
    assert any(r["sha"].startswith("835b191f") for r in steps)
    assert steps[-1]["sha"]==m["last_verified_research_commit_before_manifest"]


def test_required_independent_best_monte_carlo_and_full_frw_tasks_exist():
    m=manifest()
    priorities=m["next_required_work"]
    assert [p["priority"] for p in priorities]==[1,2,3,4]
    assert "BEST" in priorities[0]["task"]
    assert "FRW" in priorities[1]["task"]
    assert "zero-mode" in priorities[2]["task"]
    assert "CLASS" in priorities[3]["task"]


def test_analytic_number_and_energy_reference_magnitudes():
    metrics=manifest()["diagnostics"]
    assert abs(metrics["vacuum_exact_number_production_m_over_T_0p1_g1_4_Gamma_1"]-
               .0164959084596)<1e-12
    assert abs(metrics["thermal_parent_energy_loss_m_over_T_0p1_g1_4_Gamma_1"]-
               .03646818570688685)<1e-12
