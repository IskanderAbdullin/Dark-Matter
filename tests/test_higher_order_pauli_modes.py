"""Basis-size sensitivity of the Pauli thermal eigenvalue Ritz approximation."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from higher_order_pauli_modes import ritz_hierarchy
from linearized_pauli_modes import linearized_modes


def test_two_mode_basis_matches_independent_previous_galerkin_implementation():
    r=ritz_hierarchy(n=12,cutoff=24,modes=2)
    ref=linearized_modes(n=12,cutoff=24)
    np.testing.assert_allclose(r["Ritz_rates_GeV"],
                               ref["generalized_rates_GeV"],rtol=1e-10)
    assert r["number_orthogonality_max"]<1e-12
    assert r["basis_orthonormality_max"]<1e-11


def test_growing_basis_shows_thermal_single_rate_is_not_true_kinetic_gap():
    results=[ritz_hierarchy(n=20,cutoff=32,modes=k) for k in (2,3,4,5)]
    rates=[x["Ritz_rates_GeV"][0] for x in results]
    assert all(0<rates[i+1]<rates[i] for i in range(3))
    assert rates[-1]<.7*rates[0]
    assert all(r["number_orthogonality_max"]<1e-10 for r in results)
    assert all(r["basis_orthonormality_max"]<1e-10 for r in results)
    assert all(min(r["Ritz_rates_GeV"])>0 for r in results)


def test_quadrature_refined_5_mode_not_basis_converged():
    a=ritz_hierarchy(n=20,cutoff=32,modes=5)
    b=ritz_hierarchy(n=24,cutoff=32,modes=5)
    assert abs(a["Ritz_rates_GeV"][0]/b["Ritz_rates_GeV"][0]-1)<3e-4
    assert abs(a["Ritz_rates_GeV"][-1]/b["Ritz_rates_GeV"][-1]-1)<3e-4
    assert "finite Galerkin" in b["warning"]


def test_g_squared_scales_all_ritz_rates():
    a=ritz_hierarchy(n=12,cutoff=24,modes=3,G=.001)
    b=ritz_hierarchy(n=12,cutoff=24,modes=3,G=.002)
    np.testing.assert_allclose(
      np.asarray(b["Ritz_rates_GeV"])/a["Ritz_rates_GeV"],
      [4]*3,rtol=1e-12)


def test_invalid_model_inputs():
    with pytest.raises(ValueError):ritz_hierarchy(modes=6)
    with pytest.raises(ValueError):ritz_hierarchy(modes=0)
    with pytest.raises(ValueError):ritz_hierarchy(M=0.)
