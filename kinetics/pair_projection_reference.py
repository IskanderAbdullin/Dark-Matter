"""Number/energy-conserving reference projection for positive pair events."""
import numpy as np


def bracket_energy(energy, shell_mean):
    """Return adjacent indices and high-node weight; refuse extrapolation."""
    x, e = np.asarray(energy, float), np.asarray(shell_mean, float)
    if (e.ndim != 1 or e.size < 2 or not np.all(np.isfinite(e))
            or not np.all(np.diff(e) > 0) or not np.all(np.isfinite(x))):
        raise ValueError("need finite, strictly increasing shell energy means")
    if np.any(x < e[0]) or np.any(x > e[-1]):
        raise ValueError("event outside the shell-mean energy support")
    j = np.clip(np.searchsorted(e, x, side="right"), 1, e.size-1)
    w = (x-e[j-1])/(e[j]-e[j-1])
    return j-1, j, w


def deposit(energy, shell_mean):
    """Deposit one count per event; caller must separately handle Pauli limits."""
    i, j, w = bracket_energy(energy, shell_mean)
    n = len(shell_mean)
    return np.bincount(i, weights=1-w, minlength=n)+np.bincount(j, weights=w, minlength=n)
