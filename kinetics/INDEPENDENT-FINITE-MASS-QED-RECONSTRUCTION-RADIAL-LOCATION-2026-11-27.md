# Locate the bare virtual change from conservative occupation smoothing

Use the 5% requested transition reconstruction from the
[fixed shell number and energy reference](INDEPENDENT-FINITE-MASS-QED-CONSERVATIVE-VIRTUAL-RECONSTRUCTION-2026-11-27.md).
On each of the four central Born events compare its internal occupation
directly to the archived steps, holding the angular bare kernel fixed.
The external event is not at a discontinuity. Integrate the difference
rather than subtracting two independently large full-support integrals.

Split the result by absolute internal energy displacement:
0 to .0003, .0003 to .02, .02 to .25, and the remaining archived
support. In the innermost symmetric band subtract the known odd
term proportional to the central occupation difference. Its
symmetric principal value is zero. In a separate partition use
whether the uncut charged angular pole lies inside the polar aperture.
Both partitions use the exact same bare kernel and occupation difference.

The conservative-minus-step differences are −1.1673/−8.7280
(fine parent/daughter) and −12.0002/−13.3771 (third).
The innermost band contributes at most 0.18% of the net
difference in magnitude. The dominant band is .0003 to .02 for every
case, with a partly compensating contribution farther out.
About 25–78% of the signed net difference comes from outside the
polar aperture, where the uncut charged angular pole is not on the
integration interval. Thus this particular shape sensitivity is not
mainly a pointwise radial coincidence effect. Refining the direct
difference mesh changes its bands by less than 2e−10, and their
sums independently recover the original full support difference
within 1.1e−8 absolute.

This is a location audit of an isolated bare positive charged
cut, not a proof that other terms in a full hot-medium calculation
are small near the overlap. It does not fix an unknown physical line
shape or continuum occupation.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_reconstruction_radial_location_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-reconstruction-radial-location-2026-11-27.json).
