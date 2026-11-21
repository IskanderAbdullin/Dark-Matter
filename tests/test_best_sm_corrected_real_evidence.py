"""Real pinned third-party BEST+finite-radiation conservation repair evidence.

This is an archived 2026-11-21 experimental result, not a CLASS relic.
The previous experiment conserved a BGK trapezoid moment only but leaked
BEST's actual log-Simpson global number/energy. Test both old negative
result and corrected physical time-step regression.
"""
from __future__ import annotations
import json
from pathlib import Path
import math


def load(filename):
    file=Path(__file__).resolve().parents[1]/"kinetics/results"/filename
    return json.loads(file.read_text(encoding="utf-8"))


def test_both_pinned_best_runs_have_verifiable_provenance():
    old=load("best-sm-coupled-original-2026-11-21.json")
    new=load("best-sm-coupled-corrected-2026-11-21.json")
    assert old["verified_real_external_run_id"]==44925693596
    assert new["run_id"]==44929647245
    sha="c70875be714e7a6bc6d38578c34bb0dc9f5aefa9"
    assert old["pinned_BEST_git_commit"]==sha
    assert new["external_best_sha"]==sha
    assert old["negative_results"]["same_grid_BGK_vs_BEST_moments_not_identical"]
    assert new["claims"]["exact_BEST_BGK_same_discrete_measure_proven"]
    assert new["claims"]["external_BEST_AND_temperature_feedback_ACTUALLY_executed"]


def test_corrected_BGK_extremely_tight_independent_mother_number_and_heat_Q():
    new=load("best-sm-coupled-corrected-2026-11-21.json")
    for scenario in new["scenarios"]:
        assert scenario["largest_substep_number_error"]<1e-12
        assert scenario["largest_substep_dark_plus_SM_energy_Q_error"]<1e-12
        assert abs(scenario["relative_charge_phi_error"])<1e-9
        assert abs(scenario["relative_charge_chi2_error"])<1e-9
        assert abs(scenario["relative_FRW_pressure_energy_error"])<1e-7


def test_energy_pressure_work_error_regains_expected_second_order_and_avoids_plateau():
    before=load("best-sm-coupled-original-2026-11-21.json")
    after=load("best-sm-coupled-corrected-2026-11-21.json")
    rows=[r for r in after["scenarios"] if r["gamma_el"]==100000]
    assert [r["steps"] for r in rows]==[4,8,16]
    errors=[abs(r["relative_FRW_pressure_energy_error"]) for r in rows]
    assert 3.7<errors[0]/errors[1]<4.3
    assert 3.7<errors[1]/errors[2]<4.3
    assert errors[-1]<1e-9
    old_err=before["negative_results"]["gamma1e5_total_energy_error_plateau"]
    assert errors[-1]<old_err/100


def test_thermally_coupled_case_is_physical_heat_exchange_not_extra_mother_creation():
    d=load("best-sm-coupled-corrected-2026-11-21.json")
    baseline=d["scenarios"][0]
    strong=d["scenarios"][-1]
    assert baseline["gamma_el"]==0
    assert baseline["SM_comoving_energy_delta"]==0
    assert math.isclose(baseline["final_TSM"],1/math.sqrt(1.4),rel_tol=1e-14)
    assert strong["SM_comoving_energy_delta"]<0
    assert strong["final_TSM"]<baseline["final_TSM"]
    assert strong["phi_number"]>0
    assert strong["phi_f_qmin"]>1
    assert strong["phi_zero_N0"]>0
    assert not d["claims"]["microphysical_thermalization_rate_derived_from_lagrangian"]
    assert not d["claims"]["real_86eV_relic_mass_BEC_or_CLASS_validated"]
