# Preserving hard moments alone does not determine the frozen local sign

For the same fixed native packet geometry, native parent, scalar
seed, and original daughter interval as the [pointwise control](INDEPENDENT-EARLY-POLE-DAUGHTER-POINTWISE-ARCHIVE-2026-11-27.md),
allow each daughter occupation to vary between zero and one.
Within **each** retained daughter cell, fix the number and
bare energy of each arm to its own prepared moments. In the
stronger control also hold each arm's reciprocal-energy moment
at the value implied by its reconstructed occupation. Scale that
third constraint using a centered reciprocal basis rather than
subtracting nearly equal raw moments in the linear program.

The frozen paired rate is affine in the two daughter occupations.
Use a separate daughter-momentum integral for its positive
influence kernel. Split at the allowed packet boundaries, then
find piecewise-constant feasible extrema of each arm on increasingly
fine subcells. Repair moment residuals by a tiny move toward
an interior positive occupation followed by an exact numerical
moment projection. The verified feasible witnesses are not
inferred physical occupations. An independent fitted-shape
integral reproduces the prior paired packet rates within
3.23e-13. The largest checked witness moment residual is 2.23e-18.

| Hard depth × cold modes | Feasible paired range, number and energy | Feasible paired range, also reciprocal energy |
|---|---:|---:|
| 1 × 8 | [-.095001, +.079802] | [-.088439, +.061055] |
| 1 × 32 | [-.099111, +.083388] | [-.089546, +.058585] |
| 4 × 8 | [-.038788, +.017534] | [-.027664, +.005174] |
| 4 × 32 | [-.039456, +.016204] | [-.025895, +.002292] |

Thus admissible shapes compatible with the specified moments
give both signs in all four frozen configurations, including
with the extra constraint. Adjacent phases need not be smooth
in these witnesses. The result does not specify whether these
shapes are realized or whether either sign persists dynamically.
It identifies uncertainty not settled by accurate quadrature
or exact enforcement of those moments.

For an exploratory outer envelope, at fixed moment weights
(A), moments (m) and influence (W), the continuous dual
inequality for the minimum is
(lambdacdot m+intmin(W-lambdacdot A,0)).
Apply the same formula to the negated objective for the maximum.
The receipt records these integrals at two quadrature orders;
the formula is exact for exactly computed integrals, but the
reported numerical outer envelope is not an interval-arithmetic
proof. The sign conclusion uses the *inner* feasible witnesses.

Run `early_pole_daughter_frozen_moment_envelope.py` on the four
accepted hard/cold fixtures, with `--oversample 4` and 12, both
with and without `--fix-gauge`, and check with
`early_pole_daughter_frozen_moment_envelope_compare.py`.
