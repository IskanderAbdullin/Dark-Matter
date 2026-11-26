"""Run: python kinetics/test_massive_shell_moments.py > evidence.json"""
import json
import math
import platform

import numpy as np
import scipy
from scipy.integrate import quad

from massive_shell_moments import shell_weights


def rel(a, b):
    return abs(a-b) / max(abs(a), abs(b), np.finfo(float).tiny)


def reference(lo, hi, mu):
    d = hi-lo
    def moment(w):
        def integrand(t):
            q = lo + d*t
            return q*q*w(q, math.hypot(q, mu))
        return d * quad(integrand, 0, 1, epsabs=0, epsrel=3e-12)[0]
    return np.array([moment(w) for w in (
        lambda q, e: 1, lambda q, e: e,
        lambda q, e: q*q/e if e else 0, lambda q, e: mu*mu/e if e else 0)])


if __name__ == "__main__":
    rng = np.random.default_rng(251126)
    number = 1600
    quad_errors = np.zeros(4)
    identity_error = derivative_error = 0.
    cases = [(0., 1., 0.), (0., 1., 1e-8), (0., 1e-12, 1.),
             (1e-12, 1e-12*(1+1e-10), 1.), (0., 10., 1.)]
    for _ in range(number):
        q = 10**rng.uniform(-11, 4)
        cases.append((q*rng.choice([0., 1.]), q*(1 + 10**rng.uniform(-11, 2)),
                      10**rng.uniform(-12, 7)))
    h = 2e-5
    for lo, hi, mass in cases:
        edges = [lo, hi]
        n, e, p, t = vals = np.array(shell_weights(edges, mass))[:, 0]
        quad_errors = np.maximum(quad_errors, [rel(a,b) for a,b in
                                             zip(vals, reference(lo, hi, mass))])
        identity_error = max(identity_error, abs(e-p-t)/e)
        em = shell_weights(edges, mass*np.exp(-h))[1][0]
        ep = shell_weights(edges, mass*np.exp(h))[1][0]
        derivative_error = max(derivative_error, abs((ep-em)/(2*h)-t)/e)
        assert n > 0 and e > 0 and p >= 0 and t >= 0
    # Individual cells can collapse if energy differences are below machine epsilon.
    q = 1e-12
    near_rest_ratio = shell_weights([0., q], 1.)
    asymptotic_error = rel((near_rest_ratio[2]/near_rest_ratio[1])[0], 0.6*q*q)
    assert np.hypot(1., q) == 1.
    assert max(quad_errors) < 5e-9
    assert identity_error < 2e-15 and derivative_error < 1e-9
    assert asymptotic_error < 5e-15
    print(json.dumps({
        "cells": len(cases), "seed": 251126, "derivative_log_step": h,
        "max_relative_error_vs_adaptive_quad_N_R_3P_trace": quad_errors.tolist(),
        "max_abs_identity_residual_over_energy": identity_error,
        "max_abs_log_derivative_residual_over_energy": derivative_error,
        "near_rest_asymptotic_relative_error": asymptotic_error,
        "energy_edges_collapsed_for_mu_1_qmax_1e-12": True,
        "python": platform.python_version(), "numpy": np.__version__,
        "scipy": scipy.__version__}, indent=2))
