"""Stieltjes independent recurrence Ritz 2..11 mode basis-size falsification."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from stieltjes_pauli_modes import pauli_ritz
from higher_order_pauli_modes import ritz_hierarchy


def test_two_and_five_modes_match_independent_gram_schmidt_code():
    for m in (2,5):
        x=pauli_ritz(nonconstant_modes=m,n=16,cutoff=28)
        ref=ritz_hierarchy(modes=m,n=16,cutoff=28)
        np.testing.assert_allclose(x["rates_GeV"],ref["Ritz_rates_GeV"],rtol=1e-10)
        assert x["basis_orthogonality_max"]<1e-12


def test_up_to_11_modes_keep_decreasing_not_a_converged_physical_rate():
    records=[pauli_ritz(nonconstant_modes=m,n=20,cutoff=32) for m in (2,5,7,9,11)]
    rate=[r["rates_GeV"][0] for r in records]
    assert all(0<rate[i+1]<rate[i] for i in range(len(rate)-1))
    assert rate[-1]<.5*rate[0]
    assert all(r["basis_orthogonality_max"]<1e-12 for r in records)


def test_11mode_quadrature_n20_to_n24_is_not_basis_convergence():
    a=pauli_ritz(nonconstant_modes=11,n=20,cutoff=32)
    b=pauli_ritz(nonconstant_modes=11,n=24,cutoff=32)
    assert abs(a["rates_GeV"][0]/b["rates_GeV"][0]-1)<.001
    assert abs(a["rates_GeV"][-1]/b["rates_GeV"][-1]-1)<.001
    assert b["rate_means_no_spectral_gap_proof"]


def test_coupling_quadratic_and_invalid_mode_controls():
    a=pauli_ritz(nonconstant_modes=4,n=12,G=.001)
    b=pauli_ritz(nonconstant_modes=4,n=12,G=.002)
    np.testing.assert_allclose(
      np.asarray(b["rates_GeV"])/a["rates_GeV"],[4]*4,rtol=1e-11)
    with pytest.raises(ValueError):pauli_ritz(nonconstant_modes=12)
    with pytest.raises(ValueError):pauli_ritz(nonconstant_modes=0)
