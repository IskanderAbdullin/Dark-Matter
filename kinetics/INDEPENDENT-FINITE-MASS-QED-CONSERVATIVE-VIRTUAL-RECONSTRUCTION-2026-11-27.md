# Continuous occupation shapes with fixed number and energy in every cell

The [Born-weighted bare virtual test](INDEPENDENT-FINITE-MASS-QED-BORN-WEIGHTED-CHARGED-OUTER-2026-11-27.md)
compared internal archived steps to a nonconservative linear
interpolation. Replace the latter with an independently constrained
continuous reconstruction, keeping each archived shell's original
number **and** bare energy. This checks the occupation representation,
not a medium rate.

Use u=p³/3 so the shell number measure is du. Join neighboring
cells at a common occupation with cubic smoothstep transition
segments. Choose the intermediate plateau from the analytic number
constraint. On normalized cell coordinate t the added polynomial

    b(t)=t²(1−t)²(2t−1)

vanishes with its derivative at each endpoint and integrates to
zero in du. Its coefficient fixes the independent energy integral.
Subtract the center energy when evaluating both moments to avoid
large constant cancellation. Exact shell saturation fixes the
appropriate joint boundary value; if a nonsaturated cell overshoots
the Pauli bounds, decrease only its local transition width before
reevaluating both constraints. Check all piecewise polynomial
extrema in their own local coordinates and a separate dense sample.
An independent higher quadrature verifies the number and normalized
energy constraints. The largest tested errors in the normalized
cell constraints are below 1.2e−14.

At each selected external Born event integrate the whole bare
positive-frequency radial correction. With continuous internal
occupation the coincidence pole can subtract the single common
coefficient proportional to f(E_external); integrate its analytic
PV separately. Split at all reconstruction joints and both angular
aperture endpoints. The higher radial quadratures agree within
1.4e−6 absolute in the raw coefficient in these tests.

For the four untuned central cases, the 1.5–12% requested
transition widths change the bare result by roughly 0.04–3.8%
relative to the archived step result. The sign of the difference
is often the opposite of the simple nonconservative interpolation.
For example, at the fine daughter center the step result is
−606.305 and the simple linear result −533.959. The conservative
results for 1.5%, 5% and 12% requested widths are approximately
−610.299, −615.033 and −624.054. A tuned archive edge is admissible
for the continuous reconstruction without the pointwise step log;
its comparison to the step still uses the explicitly declared step
finite part. The widths of some tail cells are reduced to maintain
capacity, so a requested width is not the achieved width everywhere.

No thermal cut coefficient, damping law, physical continuum limit,
or net scalar production correction is inferred from these isolated
bare coefficients.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_conservative_virtual_reconstruction_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-conservative-virtual-reconstruction-2026-11-27.json).
