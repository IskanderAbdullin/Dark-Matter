# Restore the pure exterior interaction for a conditional Hubble clock

The [mixed live/exterior ledger](INDEPENDENT-FINITE-MASS-QED-LIVE-LOG-EXTERIOR-2026-11-27.md)
could subtract the pure exterior log energy consistently while the
scale factor was prescribed. If the dark energy density is passed
to a Hubble factor, its absolute normalization matters. Restore
the frozen exterior–exterior energy and its explicit comoving scale
work. In the same fixture, include the complementary free charged
energy below the selected hard floor and the exterior photon energy
and factorized interaction moment from the fixture's smooth photon
shape. The fermion log functional is still hard-only.

The pure exterior log depends only on eta=ln(a). Tabulate one
interpolating energy and use its derivative for the scale work;
independent direct evaluation of the energy and trace checks that
derivative. The logarithmic quadrature refinement changes those
exterior coefficients by <2.0e−9 relative in both grids. The
tabulation agrees with the direct trace within 1.5e−23 absolute.

Compare the imposed reference clock with a conditional multiplier
`(1+rho_dark/rho_visible)^−1/2` on all packet rates, using a fixed
visible comoving radiation energy pi²×106.75/30. Keep the fixed
reference exposure per unit eta, so this tests the closure and
normalization, not a cosmological calibration. The dark-to-visible
ratio is about 4.62% in the represented fixture, giving an
instantaneous clock multiplier of about 0.978. The final scalar
gain decreases by about 0.81% on both starts when applying this
clock rather than the imposed reference. All absolute energy
budgets still close after explicit scale work with positive
entropy production.

The restored pure exterior log is much smaller than the whole
dark energy. Compute its differential clock contribution directly
by a rationalized ratio rather than subtracting two almost equal
trajectories. It raises the clock by at most 1.09e−13 at sampled
states and admits a 1.14e−13 bound in this held-scale domain: the
magnitude of its comoving integral grows no faster than a² at fixed
occupations and the denominator is at least the visible energy.
Do not infer smallness of the omitted soft, time-dependent exterior
or collision model from this one small term. The test is not a
Friedmann fit or a physical dark-matter abundance.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_complete_exterior_friedmann_reference.py
```
[Raw receipt](results/finite-mass-qed-complete-exterior-friedmann-2026-11-27.json).
