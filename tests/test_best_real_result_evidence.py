"""Persist actual external BEST run provenance and do not overclaim convergence."""
from __future__ import annotations
import json
import math
from pathlib import Path


def evidence():
    f=Path(__file__).resolve().parents[1]/"kinetics"/"results"/"best-independent-2026-11-21.json"
    return json.loads(f.read_text(encoding="utf-8"))


def test_real_independent_best_run_is_pinned_and_evidence_links_remain():
    r=evidence()
    assert r["confirmed_external_run"] is True
    assert r["verified_best_sha"]=="c70875be714e7a6bc6d38578c34bb0dc9f5aefa9"
    assert r["run_id"]==44916394684
    assert r["workflow_url"].endswith("/44916394684")
    assert r["artifact_expiration_days"]==7
    assert r["mpi_ranks"]==2
    assert r["claims"]["CLASS_ready"] is False
    assert not r["claims"]["independent_seed_variance_measured"]


def test_mathematical_number_error_is_reproducible_from_actual_reported_rates():
    r=evidence()
    assert len(r["cases"])==4
    assert {x["grid_n"] for x in r["cases"]}=={28,56}
    assert {x["vegas_neval"] for x in r["cases"]}=={3000,12000}
    target=r["exact_vacuum_phi_number_rate"]
    for case in r["cases"]:
        assert math.isclose(case["phi_number_rate"]/target-1,
                            case["relative_phi_error"],rel_tol=1e-10,abs_tol=1e-12)


def test_real_2x2_data_exposes_large_coarse_grid_energy_bias():
    r=evidence()["cases"]
    lower=[x for x in r if x["grid_n"]==28]
    finer=[x for x in r if x["grid_n"]==56]
    assert all(.0039<abs(x["relative_parent_energy_error"])<.0041 for x in lower)
    assert all(abs(x["relative_parent_energy_error"])<.0001 for x in finer)
    assert all(x["relative_source_shape_rms_to_analytic_peak"]<.001 for x in r)


def test_equilibrium_net_is_small_but_not_exactly_zero():
    b=evidence()["detailed_balance"]
    assert 0<b["chi1_net_over_gross"]<1e-4
    assert 0<b["phi_net_over_gross"]<1e-4
    assert 0<b["chi2_net_over_gross"]<1e-4
    assert b["total_energy_rate"]!=0
    assert not evidence()["claims"]["closed_FRW_or_nonlinear_primordial_spectrum_validated"]
