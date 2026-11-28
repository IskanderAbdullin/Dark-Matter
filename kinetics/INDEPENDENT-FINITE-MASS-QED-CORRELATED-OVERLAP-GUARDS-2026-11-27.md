# Two-regulator paths at the bare charged/photon overlap

The [single charged regulator](INDEPENDENT-FINITE-MASS-QED-CHARGED-ABSORPTIVE-GUARD-2026-11-27.md)
leaves the photon bare. Here compare that choice with a correlated
mathematical photon invariant regulator, using the two scalar-latest
causal signs of the [contour control](INDEPENDENT-FINITE-MASS-QED-SINGLE-CUT-CONTOUR-2026-11-27.md).
No physical damping law or screening scale is used.

Let eta be the absolute charged invariant guard, g_p the uncut
charged angular denominator slope at coincidence, and A and B the
local photon angular width and charged root shift slopes. Define
d_c=eta/(|g_p| sqrt(A²+B²)) and use
x=delta_E/d_c and y=(u−c_hat)/(d_c sqrt(A²+B²)). For a fixed
dimensionless parameter lambda the absolute photon invariant guard
is lambda R_0 (A²+B²) d_c²/(1−c²), with R_0=p².
This follows quadratic photon and linear charged scaling in the overlap.

Compute the difference between two-width and charged-only amplitudes
*before* integrating. For the photon factor rationalize the difference
of the two square-root inverses. In the full chart retain the original
numerator corrections for both invariant guards. Independently compute
the local infinite-angle model: at scaled x, the photon discriminant
has the form y²+(A x/sqrt(A²+B²))² with a signed imaginary addition.
Under (x,y)→(−x,−y) its difference kernel is minus its complex conjugate
after including the scalar-latest charged denominator. The local real
shift is therefore odd and the imaginary shift even. Constant radial
occupations cancel the leading real shift, whereas a radial step keeps it.

As both invariant guards tend to zero on these correlated paths, with
lambda=.03, 1, or 30 fixed, the exact local-window chart agrees with the
scaled limiting model at the same support within 1e−8 absolute.
Extending the scaled support and independently refining both meshes
changes the limit by less than 2.4e−8 in these tests.

For lambda=1 the untuned real shifts vanish in the limit. At the
deliberately tuned fine daughter and third parent steps, the real
offsets are −4.09125 and −2.66367 in the stripped bare units; they
vary with lambda. The imaginary offset is finite also at the
untuned center, +2.19590 on the fine daughter and −2.37394 on the
third parent time-order channel. These are the differences relative to
taking the photon guard to zero first, not divergent logarithmic
coefficients or observable widths. A converged one-width limit
alone does not fix this finite subtraction at an exact archive step.

No gauge-invariant dressed self-energy, spectral overlap cancellation,
or inclusive real/virtual medium rate follows from these regulator
paths.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_correlated_overlap_guards_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-correlated-overlap-guards-2026-11-27.json).
