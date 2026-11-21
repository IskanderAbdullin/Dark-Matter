"""Regression tests for the existing benchmark's discrete collision operator."""
from __future__ import annotations
import math
import sys
from pathlib import Path

import pytest
from scipy.special import zeta

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "kinetics"))
from equilibrium_audit import condensate_criterion, discrete_equilibrium_audit


def test_condensate_requirement_is_conditional_and_exact():
    c = condensate_criterion(g_phi=1, g_chi2=4)
    assert not c["smooth_chemical_equilibrium_possible_at_C0"]
    assert math.isclose(c["formal_minimum_zero_mode_number_over_T3"],
                        2*zeta(3, 1)/math.pi**2, rel_tol=1e-13)
    # g_chi2=1 reverses the conclusion, so degeneracy assumptions matter.
    assert condensate_criterion(g_phi=1, g_chi2=1)["smooth_chemical_equilibrium_possible_at_C0"]


@pytest.mark.parametrize("mu_phi", [0.0, -1.0])
def test_continuum_balance_does_not_imply_discrete_stationarity(mu_phi):
    r = discrete_equilibrium_audit(50, 64, 48, mu_phi=mu_phi)
    assert r["original_boson_L1_rate_over_vacuum"] > 0.01
    assert r["reweighted_inverse_boson_L1_rate_over_vacuum"] < 1e-10
    assert r["reweighted_inverse_fermion_L1_rate_over_vacuum"] < 1e-10
    assert abs(r["original_daughter_number_rate_mismatch_over_vacuum"]) < 1e-11
    assert r["analytic_event_level_balance_max_abs_mu0"] < 1e-10


def test_discrete_stationarity_defect_decreases_with_resolution():
    a = discrete_equilibrium_audit(50, 64, 48)
    b = discrete_equilibrium_audit(100, 140, 72)
    assert b["original_boson_L1_rate_over_vacuum"] < a["original_boson_L1_rate_over_vacuum"]
