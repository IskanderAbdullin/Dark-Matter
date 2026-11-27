# Weight the two bare positive-frequency cuts over the Born aperture

The [external edge reference](INDEPENDENT-FINITE-MASS-QED-CHARGED-OUTER-EDGE-PROJECTION-2026-11-27.md)
used a narrow manufactured panel. Use instead the full frozen
Born support at Q=4.1e−6, with two tested endpoint exclusions.
The bare scalar Born invariant and the remaining energy measure
are constant with respect to the parent energy on this support.
This makes it possible to compare averaged virtual corrections with
the actual archived external occupation weights, without selecting
physical thermal coefficients for the internal cuts.

The stimulated forward and inverse weights are
f_P(1−f_b)(1+f_s) and f_b(1−f_P)f_s. Divide both by 1+f_s before
quadrature. Evaluate their difference as
f_P−f_b+f_b(1−f_P)/(1+f_s), avoiding cancellation of the large
scalar factor. Independently reproduce the two-way equilibrium
balance on the same kinematics to 3.4e−16 relative. The archived
f_s is about 1.07e16 on fine and 1.87e15 on third; division is
only a numerical normalization for the stated shape comparison.

Inside the local radial coincidence collar also split at any
nearby angular aperture endpoint, then stretch the logarithmic
integration panel to resolve its endpoint singularity. The
older central angular/radial reference remains unchanged when
using its default. For the new Born average, use two outer orders
on panels split at all external occupation edges. Shrink the
fractional endpoint exclusions from .002 to .0005. For the
reported averaged vertex components the largest order difference
is .0027 in raw units, and the maximum change from this exclusion
test is .0021. The profile at the exactly collinear boundary
is not assigned a pointwise value.

On the tighter support, the two archived step-reconstruction
cut averages sum to −1181.46815 (forward) and −1181.53906
(net) on fine, and −1559.51906 and −1559.35472 on third.
Changing forward-only to net weights barely changes these bare
averages in the sampled window. In contrast, changing the **internal**
occupation representation from archive steps to linear energy
interpolation changes the signed sum to −1125.69888 and
−1374.65945, about 4.7% and 11.8% in magnitude. The interpolation
is not moment preserving, and its difference is not a bound
on a continuum rate.

These coefficients are only for the two bare positive-frequency
charged cuts on the vacuum Born support. They are not the complete
thermal vertex, physical decay rate, or a justification for
multiplying a kinetic evolution by these numbers.

Replay from kinetics (performs many full radial solves):
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_born_weighted_charged_outer_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-born-weighted-charged-outer-2026-11-27.json).
