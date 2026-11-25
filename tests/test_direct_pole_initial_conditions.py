"""Independent harmonic and correlated likelihood controls."""
from pathlib import Path
import sys

import numpy as np
from scipy.special import gamma, jv

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'kinetics'))
from direct_pole_geometry import profile_ruler
from direct_pole_initial_conditions import rd_trajectory


def test_radiation_oscillator_harmonic_bessel_limit():
    _, solution = rd_trajectory(3., 0., 10., 45., rtol=1e-10)
    t = np.geomspace(.001, 45., 40)
    exact = 2**.25*gamma(1.25)*jv(.25, t)/t**.25
    assert np.max(abs(solution.sol(t)[0]-exact)) < 1e-9


def test_profile_ruler_uses_off_diagonal_covariance():
    d = np.array([2000., 3000., 4000.])
    y = d/150 + np.array([1., -1, .6])
    c = np.array([[.5, .2, .03], [.2, .8, .1], [.03, .1, .5]])
    result = profile_ruler(y, c, d)
    residual = y-d/result['ruler_Mpc']
    assert abs(d @ np.linalg.solve(c, residual)) < 1e-9
    assert np.isclose(result['chi2'], residual @ np.linalg.solve(c, residual))
