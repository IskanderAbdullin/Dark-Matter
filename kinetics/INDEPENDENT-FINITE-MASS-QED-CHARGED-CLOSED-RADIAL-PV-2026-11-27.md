# Exact positive charged-cut angles and the radial archive-edge test

This follows the [bare charged geometry](INDEPENDENT-FINITE-MASS-QED-THERMAL-CHARGED-PV-GEOMETRY-2026-11-27.md)
and its [local overlap check](INDEPENDENT-FINITE-MASS-QED-CHARGED-OVERLAP-ASYMPTOTIC-2026-11-27.md).
The positive-frequency cut is still only a bare mathematical piece in
Feynman gauge. Both charged species are kept massive.

## Exact angular primitive

Write R = pr for the product of the two charged radial momenta,
c for the external cosine, u for the polar integration variable,
and b for the nonnegative invariant mass/rapidity shift computed
without subtracting nearly equal dot products. The exact discriminant is

    X = R²(u−c)² + 2bR(1−cu) + b²
      = R²[(u−c_hat)² + a_hat²],
    c_hat = c(1+b/R),    a_hat² = (1−c²)b(2R+b)/R².

The uncut inverse photon azimuthal average is −1/(2√X). After
u−c_hat = a_hat sinh(y), the remaining polar integrals are elementary:
a difference of inverse hyperbolic sines for the plain photon average and
a difference of real logarithms for the pole integral. The latter applies
both inside the angular aperture (Cauchy principal value) and outside it
(ordinary integral). The absolute values in the logarithms are essential.
Resolved direct angular quadrature agrees within 1.3e−13 relative;
the independently regularized overlap at offsets down to 1e−8
agrees within 2.4e−12.

## Radial cutoff and archived weights

For an internal energy displacement δ the angular kernel has leading
odd term C/δ. The radial measure contributes r/2. With continuous
occupation across δ=0, the paired odd term cancels at fixed external
energy. For an external energy deliberately set to a native occupation
edge, the discontinuity leaves a radial log with coefficient
C r(0) [f(0+)−f(0−)]/2. This is a pointwise property of the archived
step representation. It is not a new in-medium divergence.

At two untuned central events per grid the log slope tends to zero.
At the tuned fine daughter edge the predicted log slope is +0.7227020
and the observed one is +0.7226865. At the tuned third-grid parent edge
the numbers are +0.6159588 and +0.6159777.
These last-decade relative defects are at most 3.1e−5.

Integrating the whole archived radial interval, we split both aperture
endpoints and all occupancy junctions and subtract the leading odd term.
Its weighted primitive is elementary on every constant or linear energy
cell. At an exact step edge the reported step number is the symmetric
finite part at reference radius 3e−4; changing that radius changes the
answer. The guard-based and analytic reconstructions agree within 9.3e−7
in the raw coefficient units, and the higher quadrature orders agree
within 2.5e−9.

For the four untuned central cuts, replacing the archived energy steps by
a continuous linear interpolation between cell centers changes the
full bare coefficient by approximately 2.8–13.5%. This is a
representation-sensitivity control, not a recommended new discrete
distribution or a conservative remap. The interpolation avoids the
pointwise step logarithm but is not guaranteed to preserve the original
moments. Neither a finite part at a sampled step nor a choice of
interpolation constitutes a physical cutoff or convergence proof.
External energy integrated against a regular weight can instead see an
integrable logarithm around an edge.

Nothing here gives the relative real-time cut weights, the matched
on-shell external residues, the corresponding real contributions, or an
inclusive medium rate. The [retarded forward-scattering construction](https://arxiv.org/abs/hep-th/9809177)
is relevant to that next step, not a justification for applying these
isolated coefficients to the cosmological collision operator.

Replay from the kinetics directory with the two existing spectrum archives:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_charged_closed_radial_pv_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-charged-closed-radial-pv-2026-11-27.json).
