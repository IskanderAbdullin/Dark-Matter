# Full-support real and imaginary charged-pole guard comparison

The [local complex guard](INDEPENDENT-FINITE-MASS-QED-CHARGED-COMPLEX-GUARD-2026-11-27.md)
checked dispersive overlap near an isolated positive cut. Here
integrate the complex angular primitive over the whole occupied
radial domain and follow both its real and its imaginary component.
The same +i epsilon convention is used in both controls. It is a
mathematical guard, not a medium spectral width.

The real log coefficient at a native occupation step is the previously
tested odd dispersive term times the jump. The imaginary angular
delta coefficient at positive or negative radial displacement has
the same leading sign for this chosen +i epsilon. If A is the
angular width slope and B the charged-root shift slope, its radial log
coefficient after weighting is

    J_im = pi (P·b) (f_+ + f_-) / [2 Q p² sqrt(A²+B²)].

Thus no occupation jump is required for an imaginary logarithm.
The untuned fine daughter and third parent absorptive coefficients are
+.128054 and +.528075, with zero real log coefficients. The tuned
coefficients are +.216708 and +.427832 in the imaginary part, with
real coefficients +.722702 and +.615959 respectively. The signs and
units concern the bare single-cut amplitudes only.

The shifted angular primitive has both dispersive and absorptive
conversion constants, obtained by replacing the real part by the
imaginary part in the corresponding normalized infinite-angle local
integral. Apply both conversions over the complete radial support
and compare independently with the sharp radial integral using
the real PV and polar imaginary delta formulas. The final slopes
are within 5.6e−7 of the analytic coefficients. Finite parts match
within 1.3e−6 absolute, with the two radial meshes agreeing to 2.1e−9
or better. The comparison splits occupation edges and aperture
endpoints and includes geometrically refined exterior panels.

This result does not sum all causal cuts. In particular the two
physical charged-leg conventions need not use the same angular
imaginary sign. A bare single-cut imaginary contribution is not a
measurable damping rate and does not determine a photon regulator.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_charged_absorptive_guard_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-charged-absorptive-guard-2026-11-27.json).
