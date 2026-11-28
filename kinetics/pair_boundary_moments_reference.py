"""Keep number and energy of edge events as explicit residual moments.

This is diagnostic accounting, not an evolved distribution. Boundary moments
are not FV occupations and do not define Pauli blocking or inverse rates.
"""
import numpy as np

from pair_projection_reference import bracket_energy


def project_with_edge_moments(energy, weights, means, energy_edges):
    """Return interior counts and [low, high] (count, energy) moment pairs."""
    x, w, m, edges = map(lambda a: np.asarray(a, dtype=float),
                         (energy, weights, means, energy_edges))
    if (x.ndim != 1 or w.shape != x.shape or m.ndim != 1 or len(m) < 2
            or edges.shape != (2,) or not all(np.all(np.isfinite(a))
                                                for a in (x, w, m, edges))
            or not np.all(np.diff(m) > 0) or not edges[0] < m[0]
            or not m[-1] < edges[1]):
        raise ValueError("invalid event arrays or energy grid")
    if np.any((x < edges[0]) | (x > edges[1])):
        raise ValueError("event outside the physical grid; extend the grid")
    low, high = x < m[0], x > m[-1]
    ok = ~(low | high)
    i, j, z = bracket_energy(x[ok], m)
    count = (np.bincount(i, weights=w[ok]*(1-z), minlength=len(m))
             + np.bincount(j, weights=w[ok]*z, minlength=len(m)))
    residual = np.array([(np.sum(w[mask]), np.sum(w[mask]*x[mask]))
                         for mask in (low, high)])
    return count, residual
