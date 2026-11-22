"""Logit-interpolated Pauli events repair bin entropy algebra, not accuracy."""
import math
import numpy as np
import pytest
from scipy.special import expit
from kinetics.paired_chi_lepton_logit_galerkin_5d import (
    logit_interpolant,entropy_consistent_projection,
)
from kinetics.nonthermal_chi_lepton_5d_weak import weak_moments


def oscillatory():
    return (lambda E:expit(-np.asarray(E)+.3*np.sin(4*np.asarray(E))),
            lambda E:expit(-np.asarray(E)+.3*np.cos(4*np.asarray(E))))


@pytest.mark.parametrize("bins",(10,40,160))
def test_same_logit_interpolation_and_deposition_gives_exact_positive_entropy(bins):
    f,g=oscillatory()
    x=entropy_consistent_projection(f,g,bins=bins)
    assert x["direct_event_entropy"]>0
    assert x["interpolated_bin_entropy"]>0
    assert abs(x["entropy_projection_relative_error"])<5e-11
    assert abs(x["number_chi_residual"])<2e-24
    assert abs(x["number_lepton_residual"])<2e-24
    assert math.isclose(x["weak_Qchi"],x["projected_Qchi"],rel_tol=5e-12)
    assert math.isclose(x["weak_Qchi"],-x["projected_Qlepton"],rel_tol=5e-12)


def test_entropy_safety_does_not_remove_severe_coarse_nonlinear_shape_aliasing():
    f,g=oscillatory()
    continuum=weak_moments(f,g)
    low=entropy_consistent_projection(f,g,bins=10)
    fine=entropy_consistent_projection(f,g,bins=320)
    assert low["direct_event_entropy"]/continuum["entropy_production"]<.1
    assert abs(fine["direct_event_entropy"]/continuum["entropy_production"]-1)<.02
    assert low["direct_event_entropy"]>0


def test_affine_fd_log_odds_exactly_reproduced_between_energy_bins():
    grid=np.linspace(.1,24.1,21)
    f=lambda E:expit(.12-np.asarray(E)/.9)
    interpolated=logit_interpolant(grid,f(grid))
    energies=np.linspace(.1,24.1,41)
    assert np.max(abs(interpolated(energies)-f(energies)))<2e-14
    l=lambda E:expit(-.08-np.asarray(E))
    q=entropy_consistent_projection(f,l,bins=20)
    direct=weak_moments(f,l)
    assert math.isclose(q["weak_Qchi"],direct["Qchi"],rel_tol=2e-11)
    assert math.isclose(q["direct_event_entropy"],
                        direct["entropy_production"],rel_tol=2e-11)


def test_reject_nonphysical_grid_or_fermi_endpoint():
    with pytest.raises(ValueError):
        logit_interpolant(np.array([0,1,2,3]),np.array([.1,.2,1.,.3]))
    with pytest.raises(ValueError):
        logit_interpolant(np.array([0,2,1,3]),np.array([.1,.2,.3,.4]))
