"""Actual independently run BEST massive scalar source evidence, not fake BEC."""
from __future__ import annotations
import math
import json
from pathlib import Path


def data():
    return json.loads((Path(__file__).resolve().parents[1]/
        "kinetics/results/best-massive-phi-real-2026-11-21.json").read_text(encoding="utf-8"))


def test_pinned_actual_best_massive_phi_mpi_evidence_and_no_class_claim():
    d=data()
    assert d["external_best_actual_run_id"]==44917857964
    assert d["verified_true_mpi_execution"]
    assert d["git_head_in_external_best"]=="c70875be714e7a6bc6d38578c34bb0dc9f5aefa9"
    assert d["mpi_ranks"]==2
    assert d["finite_qmin_is_not_zero_mode"]
    assert d["external_BEST_initial_massive_daughter_source_confirmed"]
    assert not d["BEC_condensate_dynamics_validated"]
    assert not d["CLASS_relic_86eV_validated"]


def test_actual_independent_two_mass_two_grid_scan_is_reproducible():
    d=data()
    assert len(d["results"])==4
    assert {x["mphi"] for x in d["results"]}=={.001,.010}
    assert {x["n"] for x in d["results"]}=={56,96}
    reference=d["analytic_daughter_number_rate"]
    for x in d["results"]:
        assert math.isclose(x["phi_rate"]/reference-1,x["rel_num"],
                            rel_tol=1e-9,abs_tol=1e-13)
        assert abs(x["rel_num"])<.0001
        assert abs(x["energy_sum"])<2e-6
        assert x["shape_rms_peak"]<.0001
        assert x["best_f_dot_qmin"]>0
        assert x["analytic_f_dot_zero"]>0


def test_real_best_qmin_source_matches_exact_finite_mass_rate():
    for mass in (.001,.010):
        group=[r for r in data()["results"] if r["mphi"]==mass]
        refined=next(r for r in group if r["n"]==96)
        coarse=next(r for r in group if r["n"]==56)
        qmin_relative_error=(
            refined["best_f_dot_qmin"]/refined["analytic_f_dot_qmin"]-1)
        assert abs(qmin_relative_error)<1e-4
        assert abs(refined["rel_num"])<abs(coarse["rel_num"])
        assert abs(refined["shape_rms_peak"])<abs(coarse["shape_rms_peak"])
        assert refined["analytic_f_dot_zero"]!=refined["analytic_f_dot_qmin"]


def test_energy_flux_precision_improves_with_resolution_without_claim_of_full_energy():
    r=data()["results"]
    for mass in (.001,.010):
        base=next(x for x in r if x["mphi"]==mass and x["n"]==56)
        fine=next(x for x in r if x["mphi"]==mass and x["n"]==96)
        assert abs(fine["energy_sum"])<abs(base["energy_sum"])
        assert abs(fine["rel_phi_energy"])<abs(base["rel_phi_energy"])
