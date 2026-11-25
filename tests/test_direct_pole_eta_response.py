"""Independent angular and dynamical cold limits for eta polarization."""
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import zeta

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'kinetics'))
from direct_pole_eta_stability import EtaResponse, angular_response


def test_angular_response_against_direct_integral():
    for y in (1e-8, .002, .4, 10., 1e3):
        value = quad(lambda u: (y*u)**2/(1+(y*u)**2), 0., 1., epsabs=1e-24)[0]
        assert np.isclose(angular_response(y), value, rtol=5e-11, atol=1e-26)


def test_eta_cold_response_against_pressureless_limit():
    mass, temp, g, k, growth = 1., 1e-5, .2, .3, .1
    eta = EtaResponse(26, 1e-12)
    _, v, weights = eta.coefficients(mass, temp, g, 1)
    measured = weights @ angular_response(k*v/growth)
    number = temp**3*zeta(3)/np.pi**2
    expected = g*g*number/mass*(k/growth)**2
    assert np.isclose(measured, expected, rtol=1e-7, atol=0.)
