"""Actual RD-FRW BEST three-cutoff nonlinear massive Bose early-time no-BEC result."""
from __future__ import annotations
import json
import math
from pathlib import Path


def results():
    p=Path(__file__).resolve().parents[1]/"kinetics/results/best-rd-frw-qmin-real-2026-11-21.json"
    return json.loads(p.read_text(encoding="utf-8"))


def test_actual_pinned_external_full_frw_bose_kinetic_run():
    r=results()
    assert r["actual_external_BEST_run"]==44924098843
    assert r["verified_external_BEST_git_sha"]=="c70875be714e7a6bc6d38578c34bb0dc9f5aefa9"
    assert r["raw_artifact_id"]==14714772886
    assert len(r["cases"])==3
    assert r["numerically_converged_qmin_at_fixed_grid_and_timestep"]
    assert not r["macroscopic_BEC_at_tGamma_4e5_detected"]
    assert not r["late_time_BEC_absence_proven"]
    assert not r["SM_thermalization_or_86eV_cosmological_relic_CLASS_validated"]


def test_massive_phi_subgrid_zero_mode_has_regular_qmin_cubed_volume_scaling():
    r=results()
    qcut=[x["qmin"] for x in r["cases"]]
    assert qcut==[1e-4,5e-5,2e-5]
    m=[x["measured_N0"] for x in r["cases"]]
    assert math.isclose(m[0]/m[1],8.,rel_tol=5e-5)
    assert math.isclose(m[1]/m[2],15.625,rel_tol=5e-5)
    assert r["scaling"]["relative_N0_over_qmin_cubed_scatter"]<5e-5
    assert m[0]>m[1]>m[2]>0


def test_full_phi_number_stable_under_fivefold_cutoff_change():
    r=results()
    assert r["scaling"]["phi_number_relative_range"]<2e-6
    for c in r["cases"]:
        assert math.isclose(c["measured_N0"]/c["qmin"]**3,
                            c["N0_divided_by_qmin_cubed"],rel_tol=1e-13)
        assert abs(c["relative_chi1_phi_charge_residual"])<1e-8
        assert abs(c["relative_FRW_energy_work_residual"])<1e-8


def test_analytical_regular_cell_number_is_close_but_not_equal():
    r=results()
    for c in r["cases"]:
        assert math.isclose(
            c["measured_N0"]/c["analytic_regular_below_qmin"],
            c["measured_to_analytic_N0_ratio"],rel_tol=1e-12)
        assert abs(c["measured_to_analytic_N0_ratio"]-1)<.01
        assert abs(c["phi_qmin_f"]/c["analytic_phi_qmin_f"]-1)<.006
    assert not r["full_ir_time_Nq_Vegas_convergence_proven"]
