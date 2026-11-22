"""Eventwise conservative energy projection is NOT entropy safe on coarse bins."""
import math
import numpy as np
from scipy.special import expit
from kinetics.paired_chi_lepton_barycentric_5d import projected_pair
from kinetics.nonthermal_chi_lepton_5d_weak import weak_moments


def _smooth():
    return (lambda E:expit(-np.asarray(E)/.9+.15*np.exp(-np.asarray(E)/.9)),
            lambda E:expit(-np.asarray(E)+.05*np.exp(-np.asarray(E))))


def _oscillatory():
    return (lambda E:expit(-np.asarray(E)+.3*np.sin(4*np.asarray(E))),
            lambda E:expit(-np.asarray(E)+.3*np.cos(4*np.asarray(E))))


def test_weak_vs_binned_rate_exact_number_energy_for_massless_and_massive_lepton():
    for m in (0.,.5):
        f,g=_smooth()
        x=projected_pair(f,g,m=m,bins=20)
        weak=weak_moments(f,g,m=m)
        assert math.isclose(x["weak_Qchi"],weak["Qchi"],rel_tol=2e-11)
        assert math.isclose(x["projected_Qchi"],weak["Qchi"],rel_tol=2e-11)
        assert math.isclose(x["projected_Qlepton"],-weak["Qchi"],rel_tol=2e-11)
        assert abs(x["number_chi_residual"])<1e-24
        assert abs(x["number_lepton_residual"])<1e-24
        assert math.isclose(x["direct_event_entropy"],weak["entropy_production"],
                            rel_tol=2e-11)


def test_fd_affine_logodds_entropy_projection_matches_direct_event_exactly():
    f=lambda E:expit(-np.asarray(E)/.9)
    g=lambda E:expit(-np.asarray(E))
    for m in (0.,.5):
        x=projected_pair(f,g,m=m,bins=10)
        assert math.isclose(x["interpolated_bin_entropy"],
                            x["direct_event_entropy"],rel_tol=3e-11)


def test_conservative_projection_can_generate_negative_H_theorem_on_coarse_grid():
    f,g=_oscillatory()
    x=projected_pair(f,g,bins=10)
    assert x["direct_event_entropy"]>6e-13
    assert x["interpolated_bin_entropy"]<0
    assert math.isclose(x["weak_Qchi"],x["projected_Qchi"],rel_tol=5e-12)
    assert abs(x["number_chi_residual"])<1e-24
    assert abs(x["number_lepton_residual"])<1e-24
    finer=projected_pair(f,g,bins=320)
    assert finer["interpolated_bin_entropy"]>0
    assert abs(finer["entropy_projection_relative_error"])<.02
