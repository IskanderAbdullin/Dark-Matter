"""Comoving moments for a piecewise-constant occupation on momentum shells.

All momenta and mass are in units of the same Tref proportional to 1/a.
Weights omit the common species degeneracy factor g/(2*pi**2).
No quadrature across shell boundaries or division by near-rest energy intervals.
"""
import numpy as np

_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(32)


def shell_weights(q_edges, mu):
    """Return number, energy, three times pressure, and trace weights."""
    q = np.asarray(q_edges, dtype=float)
    if (q.ndim != 1 or q.size < 2 or not np.all(np.isfinite(q))
            or q[0] < 0 or not np.all(np.diff(q) > 0)
            or not np.isfinite(mu) or mu < 0):
        raise ValueError("require finite increasing nonnegative edges and mass")
    a, b = q[:-1], q[1:]
    d = b - a
    number = d * (b*b + a*b + a*a) / 3
    if mu == 0:
        energy = d * (a+b) * (a*a+b*b) / 4
        return number, energy, energy.copy(), np.zeros_like(energy)
    # Subdivide broad cells across the relativistic transition; no independent
    # subdivisions of the pressure and energy moments (keep a common rule).
    energy, pressure3, trace = (np.zeros_like(a) for _ in range(3))
    left = a.copy()
    live = np.arange(a.size)
    while live.size:
        lo = left[live]
        hi = b[live]
        cut = np.minimum(hi, 8 * np.maximum(lo, mu))
        half = (cut - lo) / 2
        p = lo[:, None] + half[:, None] * (_NODES[None, :] + 1)
        eps = np.hypot(p, mu)
        m = half[:, None] * _WEIGHTS[None, :] * p**2
        energy[live] += np.sum(m * eps, axis=1)
        pressure3[live] += np.sum(m * (p/eps) * p, axis=1)
        trace[live] += np.sum(m * (mu/eps) * mu, axis=1)
        left[live] = cut
        live = live[cut < hi]
    return number, energy, pressure3, trace
