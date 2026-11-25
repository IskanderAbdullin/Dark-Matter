"""Check physical normalization and stable angular integration of the NR pole."""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'kinetics'))
from direct_pole_kinetics import event_geometry, pole_average
from elastic_implicit_trial import ImplicitElastic


def test_pole_average_matches_direct_quadrature():
    for low, high, total in [(0.2, 2.0, 4.1), (5., 5.00000001, 26.), (0., 10., 100.)]:
        for b in (0.01, 7., 1e5, 1e18):
            ref = quad(lambda u: (b/(b+total-(low+(high-low)*u)**2))**2,
                       0., 1., epsabs=1e-14, epsrel=1e-11)[0]
            assert np.isclose(pole_average(total, low, high, b), ref, rtol=2e-11)
    assert np.isclose(pole_average(100., 0., 10., 1e18), 1., rtol=0., atol=1e-15)


def test_dilute_pair_normalization_against_maxwell_relative_speed():
    model = ImplicitElastic(np.r_[0., np.linspace(.1, 9., 90)], 1e10)
    f = np.exp(-model.p**2/2)
    rate = model.rate*(2*model.m)**3
    number = model.vol @ f
    pair_norm = rate*f[model.i]*f[model.j]/(32*np.pi**4*number**2)
    for b in (10., 1e10):
        measured = pair_norm @ pole_average(*event_geometry(model), b)
        # Independent relative speed for two Gaussian vectors with unit
        # one-dimensional dispersion: P(w)=w^2 exp(-w^2/4)/(2 sqrt(pi)).
        ref = quad(lambda w: w**3*np.exp(-w*w/4)/(2*np.sqrt(np.pi))*(b/(b+w*w))**2,
                   0., np.inf, epsrel=1e-10)[0]
        assert np.isclose(measured, ref, rtol=.005)
