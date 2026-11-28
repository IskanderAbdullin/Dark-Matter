# Invariant-energy weights behind the two-bath differential slope crossings

At the common FD state, partition the original positive quadrature by
invariant pair energy around the vector pole. Use exact cut points
at both ±30 and ±100 Breit–Wigner widths. Below, near, and above
contributions independently add back to the original energy moments
and the positive blocking derivative weights. The window is a diagnostic,
not a split into different reactions.

For ±100 widths, the in-window share of the number weight M0 is only
0.04236% at visible T=.30 GeV, and the corresponding positive
reciprocity weight S0 share is 3.68e−8 as a fraction. At .64135 GeV,
where the measured differential order has changed, these shares
are 98.3212% and 85.1036%. For ±30 widths the corresponding latter
figures remain 97.5524% and 84.3972%. The slopes weight the same
lab energy by different measures; strong invariant pole dominance
of M0 does not ensure the same pole fraction in the Pauli derivative.
Both quadrature orders reproduce all reported fractions.

This is a numerical explanation at the same stationary samples,
not a decomposition of a physical freeze-out evolution.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_reciprocity_spectral_split_reference.py
```
