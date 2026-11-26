"""Actual external BEST non-equilibrium Bose/Pauli evidence, fail-closed science."""
from __future__ import annotations
import json
import math
from pathlib import Path


def data():
    return json.loads((Path(__file__).resolve().parents[1]/
                       "kinetics/results/best-affinity-real-2026-11-21.json").read_text(
                       encoding="utf-8"))


def test_real_external_non_equilibrium_source_is_pinned_and_not_claimed_as_frw():
    d=data()
    assert d["verified_run"]==44916922599
    assert d["run_url"].endswith("44916922599")
    assert d["actual_best_sha"]=="c70875be714e7a6bc6d38578c34bb0dc9f5aefa9"
    assert d["MPI_ranks"]==2
    assert d["scientific_convergence_confirmed"] is False
    assert d["CLASS_ready"] is False
    assert not d["conditions"]["full_kinetic_time_evolution"]


def test_four_field_analytic_energy_moments_and_source_rates_are_consistent():
    d=data()
    ref=d["analytical_reference_default"]
    assert math.isclose(ref["number_rate"],.0008565027253502168,rel_tol=1e-9)
    assert abs(sum(ref["energy_rates"].values()))<1e-12
    assert math.isclose(ref["entropy_production"],.05*ref["number_rate"],rel_tol=1e-12)
    for r in d["cases"]:
        assert math.isclose(r["best_nphi"]/r["analytic_nphi"]-1,
                            r["relative_nphi_error"],rel_tol=1e-8,abs_tol=1e-12)
        assert math.isclose(
            r["best_energy_chi1"]+r["best_energy_phi"]+r["best_energy_chi2"],
            r["physical_energy_residual"],rel_tol=1e-8,abs_tol=1e-12)
        assert r["best_entropy_rate"]>0


def test_higher_q_resolution_reduces_number_and_energy_errors_but_not_to_zero():
    cases=data()["cases"]
    q56=next(x for x in cases if x["n"]==56 and x["mu1"]==-.15)
    q96=next(x for x in cases if x["n"]==96)
    assert abs(q96["relative_nphi_error"])<abs(q56["relative_nphi_error"])
    assert abs(q96["physical_energy_residual"])<abs(q56["physical_energy_residual"])
    assert q96["shape_rms_over_peak"]<q56["shape_rms_over_peak"]
    # The entropy estimate is more sensitive to nonzero numerical
    # physical energy/charge residuals than the net particle source.
    relative_entropy=abs(q96["best_entropy_rate"]/
                         q96["analytic_entropy_rate"]-1)
    assert .005<relative_entropy<.007
    assert abs(q96["relative_nphi_error"])<.0001


def test_reverse_nonzero_bose_channel_changes_sign_without_negative_entropy():
    r=next(x for x in data()["cases"] if x["mu1"]==-.25)
    assert r["best_nphi"]<0
    assert r["analytic_nphi"]<0
    assert r["best_nchi1"]>0
    assert r["best_nchi2"]<0
    assert r["best_entropy_rate"]>0
    assert r["analytic_entropy_rate"]>0
