"""Independent spinor trace and retarded normalization of the free charged cut."""

import gzip
import json
from pathlib import Path

import numpy as np

from early_pole_hard_yukawa_photoproduction_reference import norm_spinors, slash
from finite_mass_qed_frozen_transverse_cut_reference import response
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import PREF, Species
from finite_mass_qed_spacelike_exact_statistical_aperture_reference import cell, weight


def run():
    trace_defect = 0.0
    retarded = []
    for grid in ("fine", "third"):
        root = Path(__file__).parent / "results"
        with gzip.open(root / f"early-pole-spectrum-{grid}-m06-2026-11-26.json.gz") as stream:
            data = json.load(stream)
        parent = data["x_end"]
        charges = [Species(m, data[f"endpoint_{n}_momentum_edges"], data[f"endpoint_{n}_occupations"])
                   for m, n in ((parent, "parent"), (parent * data["mass_ratio"], "fermion"))]
        discrepancy = 0.0
        for k in (.0002, .0008, .00256, .01):
            for u in (.0047, .25, .9, .99):
                w = k * u
                s = (k - w) * (k + w)
                D = sum(cell(q.m, k, w, q.e, q.f)[2] for q in charges)
                exact = PREF * np.pi * D / (16 * k)
                soft = response(charges, u)[1]
                if min(exact, soft) <= 0:
                    raise ArithmeticError("Nonpositive retarded cut")
                discrepancy = max(discrepancy, abs(exact / soft - 1))
                for q in charges:
                    Emin, weight_at = weight(q.m, w, k)
                    E = Emin + .1
                    pz = (w * E - s / 2) / k
                    p_perp = np.sqrt((E-q.m)*(E+q.m)-pz*pz)
                    p = np.array([E, p_perp, 0., pz])
                    p_next = p + np.array([w, 0., 0., k])
                    x, y = slash(np.array([0., 1., 0., 0.])), slash(np.array([0., 0., 1., 0.]))
                    value = .5 * sum(norm_spinors(p, p_next, q.m, q.m, e) for e in (x, y))
                    trace_defect = max(trace_defect, float(abs(value / weight_at(E) - 1.)))
        retarded.append({"grid": grid, "exact_vs_soft_retarded_cut_max_rel": float(discrepancy)})
    if trace_defect > 5e-13 or max(r["exact_vs_soft_retarded_cut_max_rel"] for r in retarded) > 1.2e-5:
        raise ArithmeticError((trace_defect, retarded))
    return {"independent_spinor_trace_max_rel": trace_defect, "retarded": retarded}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
