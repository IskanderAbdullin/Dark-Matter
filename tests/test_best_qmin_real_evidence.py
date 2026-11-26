"""Independent pinned BEST non-condensate short-time zero-mode qmin^3 evidence."""
from __future__ import annotations
import json
import math
from pathlib import Path


def observed():
    return json.loads((Path(__file__).resolve().parents[1]/
       "kinetics/results/best-qmin-real-2026-11-21.json").read_text(encoding="utf-8"))


def test_actual_external_best_cutoff_scan_and_no_late_bec_claim():
    d=observed()
    assert d["real_external_best_action_run"]==44920000218
    assert d["verified_best_git_sha"]=="c70875be714e7a6bc6d38578c34bb0dc9f5aefa9"
    assert d["macroscopic_zero_mode_BEC_detected_at_this_early_time"] is False
    assert d["absence_of_future_late_time_condensate_proven"] is False
    assert d["primordial_86eV_dark_matter_or_CLASS_validated"] is False


def test_measured_massive_zero_mode_scales_as_qmin_cubed_to_0p02_percent():
    d=observed();rows=d["cases"]
    assert len(rows)==3
    assert [x["qmin"] for x in rows]==[1e-4,5e-5,2e-5]
    a,b,c=rows
    assert math.isclose(a["N0"]/b["N0"],8,rel_tol=.0003)
    assert math.isclose(b["N0"]/c["N0"],15.625,rel_tol=.0003)
    assert d["volume_ratio_checks"]["N0_divided_by_qmin_cubed_relative_scatter"]<.0003
    assert a["N0"]>b["N0"]>c["N0"]>0


def test_finite_q0_occupation_stabilizes_while_subcell_number_vanishes():
    d=observed()
    phi=[x["phi_number"] for x in d["cases"]]
    assert max(phi)/min(phi)-1<.0001
    assert all(1.9<x["phi_f_qmin"]<2.0 for x in d["cases"])
    assert all(abs(x["charge_rel"])<1e-8 for x in d["cases"])
    assert all(abs(x["energy_rel"])<1e-9 for x in d["cases"])
    assert d["cases"][-1]["N0"]<3e-16
