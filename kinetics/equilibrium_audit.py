"""Audit discrete 1<->2 detailed balance and the fixed-bath condensation threshold.

This is a diagnostic of the existing benchmark, NOT a replacement kinetic model.
No integration to finite time or cosmological claims are made here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import expit, zeta

from one_to_two_fixed_bath import FlatThermalDecay


def condensate_criterion(g_phi: int = 1, g_chi2: int = 4) -> dict:
    """Formal T=1, m=0 equilibrium with mu_parent=0 and n_phi-n_chi2=0.

    mu_chi2=-mu_phi and mu_phi<=0. For mu_phi<0, the smooth-boson
    density decreases while the fermion density increases. Consequently,
    if at mu_phi=0 n_chi2>n_phi, smooth equilibria cannot satisfy C=0.
    The surplus is the *required zero-mode occupation* IF equilibrium with
    a condensate is attainable; it is not a prediction of its formation.
    """
    z3 = float(zeta(3.0, 1.0))
    nb = g_phi * z3 / math.pi**2
    nf = (3.0/4.0) * g_chi2 * z3 / math.pi**2
    return {
        "boson_degeneracy": g_phi,
        "fermion_degeneracy": g_chi2,
        "thermal_boson_number_over_T3_at_mu0": nb,
        "thermal_fermion_number_over_T3_at_mu0": nf,
        "critical_conserved_difference_C_over_T3": nb-nf,
        "initial_difference_C_over_T3": 0.0,
        "formal_minimum_zero_mode_number_over_T3": max(0.0, nf-nb),
        "smooth_chemical_equilibrium_possible_at_C0": nf <= nb,
        "caveat": "Only a fixed-T massless ideal quantum gas and conserved daughter difference; no dynamical condensate formation is demonstrated.",
    }


def discrete_equilibrium_audit(nq: int, np1: int, nmu: int, mu_phi: float = 0.0) -> dict:
    """Evaluate the actual cell ODE on an analytic quantum equilibrium.

    Event-level balance is not sufficient: four daughter-bin pairs
    have shell-mean energy sums generally different from parent E1.
    A reweighted inverse term, using exp(e_phi_bin+e_chi_bin), is included
    strictly as a zero-residual diagnostic, not validated physics.
    """
    k = FlatThermalDecay(nq=nq, np1=np1, nmu=nmu)
    fp = 1.0 / np.expm1(k.avg_q - mu_phi)
    fc = expit(-mu_phi - k.avg_q)  # mu_chi2 = -mu_phi
    y = np.r_[fp, fc]
    dy = k.rhs(0.0, y, stimulated=True)
    rate_phi = dy[:nq] * k.phi_volume
    rate_chi = dy[nq:] * k.chi_volume
    vacuum_rate = float(np.dot(k.rateprefactor, k.parent_fermi))

    # Keep forward rates unchanged and only alter the inverse coefficient.
    # This is a diagnostic: it changes the finite-grid physical kernel.
    adjusted_phi = np.zeros(nq)
    adjusted_chi = np.zeros(nq)
    for i in (0, 1):
        ip = k.ip + i
        wp = k.wp if i else (1.0 - k.wp)
        b = fp[ip]
        for j in (0, 1):
            ic = k.ic + j
            wc = k.wc if j else (1.0 - k.wc)
            f = fc[ic]
            a = k.parent_fermi
            stat = a * ((1.0 - f) * (1.0 + b) - np.exp(k.avg_q[ip] + k.avg_q[ic]) * f * b)
            events = k.rateprefactor * wp * wc * stat
            adjusted_phi += np.bincount(ip, weights=events, minlength=nq)
            adjusted_chi += np.bincount(ic, weights=events, minlength=nq)

    event_balance, event_energy = k.detailed_balance_error()
    return {
        "grid": [nq, np1, nmu],
        "mu_phi_over_T": mu_phi,
        "mother_vacuum_forward_rate_over_T3_per_Gamma": vacuum_rate,
        "original_net_number_rate_over_vacuum": float(rate_phi.sum()/vacuum_rate),
        "original_boson_L1_rate_over_vacuum": float(np.abs(rate_phi).sum()/vacuum_rate),
        "original_fermion_L1_rate_over_vacuum": float(np.abs(rate_chi).sum()/vacuum_rate),
        "original_daughter_number_rate_mismatch_over_vacuum": float((rate_phi.sum()-rate_chi.sum())/vacuum_rate),
        "reweighted_inverse_boson_L1_rate_over_vacuum": float(np.abs(adjusted_phi).sum()/vacuum_rate),
        "reweighted_inverse_fermion_L1_rate_over_vacuum": float(np.abs(adjusted_chi).sum()/vacuum_rate),
        "analytic_event_level_balance_max_abs_mu0": float(event_balance),
        "event_energy_error_max_abs": float(event_energy),
    }


def run(quick: bool = False) -> dict:
    grids = [(50, 64, 48)] if quick else [
        (50, 64, 48), (100, 140, 72), (150, 180, 96), (180, 240, 120)
    ]
    cases = [discrete_equilibrium_audit(*grid) for grid in grids]
    for row in cases:
        assert abs(row["original_daughter_number_rate_mismatch_over_vacuum"]) < 1e-11, row
        assert row["reweighted_inverse_boson_L1_rate_over_vacuum"] < 1e-10, row
        assert row["reweighted_inverse_fermion_L1_rate_over_vacuum"] < 1e-10, row
        assert row["analytic_event_level_balance_max_abs_mu0"] < 1e-10, row
    return {
        "baseline_commit": "ea6ec461997600d5f00bf244ea808607d40ef8dd",
        "scope": "Minkowski fixed mother bath, 1<->2, massless daughters, 2026-11-21",
        "thermodynamic_criterion": condensate_criterion(),
        "discrete_operator_audit_at_mu0": cases,
        "interpretation": "The continuum event balance vanishes but the bin-evaluated ODE does not annihilate BE/FD equilibrium. Its defect shrinks on finer grids. An inverse-term reweighting exactly annihilates the thermal state, but is not a validated physical solver. A separate daughter-number invariant makes a regular massless equilibrium impossible at C=0 for g_phi=1,g_chi2=4; zero-mode/IR treatment requires dedicated kinetic development.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--output", default="outputs/kinetic_equilibrium_audit.json")
    args = parser.parse_args()
    answer = run(quick=args.quick)
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(answer, indent=2), encoding="utf-8")
    print(json.dumps(answer, indent=2))
