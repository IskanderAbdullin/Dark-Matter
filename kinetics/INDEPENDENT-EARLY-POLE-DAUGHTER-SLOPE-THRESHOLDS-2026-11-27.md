# Separate the slope required by the stored moments from the slope required for a positive witness

The [continuous aperture test](INDEPENDENT-EARLY-POLE-DAUGHTER-CONTINUOUS-SIGN-WITNESS-2026-11-27.md)
ties adjacent nodal occupancies but places no upper limit on
their slope. Here use continuous piecewise-linear occupations
on a uniform subgrid within each hard daughter cell. Keep
cell number and bare energy, and optionally the reconstructed
reciprocal-energy moment, exactly within numerical tolerances.
Set (|f_{i+1}-f_i|le L(p_{i+1}-p_i)), where momentum is
expressed in p/T. Repair a finite candidate using an interior
anchor and explicitly recheck occupations, moments and slope.
The source remains the same frozen initial paired observable.

Minimize the shared cap (L) first for feasibility of the two
sets of moments, then separately when requiring a positive
rate margin of (10^{-5}) in the finite optimization.
The two problems are not interchangeable. Evaluate the
weak dual from the raw matrices and signed optimization
multipliers, rather than relying on solver status alone.
Adjacent cell averages and their mean momenta also give an
analytic lower floor on feasibility, which does not require
a finite profile. The verified finite witness values below
come from 128 subdivisions per hard cell:

| Hard depth × cold modes | Feasible slope, three moments | Slope of positive witness, three moments | Slope of positive witness, two moments |
|---|---:|---:|---:|
| 1 × 8 | 7.032 | 11.059 | 4.230 |
| 1 × 32 | 7.032 | 10.874 | 4.230 |
| 4 × 8 | 28.215 | 624.646 | 52.532 |
| 4 × 32 | 28.215 | 598.070 | 50.792 |

If the two profile endpoints also match the fixed occupations
immediately outside the interval, the depth-four feasibility
slope rises to 46.081. The positive-witness slopes at that
depth change by less than 0.002. These are numerical finite
interpolant optima or verified near-optimal examples, not
established physical bounds on the slope of the daughter
distribution. The fitted maximum slope measured *within* a
hard cell is not the slope of a continuous reconstruction:
the fitted profiles can jump across cell boundaries.

For a separate approximate outer check, interpolate any
Lipschitz occupation by its chords. In an interval of width
(h), at fractional position (t), the absolute error is
bounded by (2Lh,t(1-t)) by the triangle inequality.
Integrate this bound against the absolute moment/source
weights. Relax the finite moment constraints by the resulting
errors and add the allowed source error to the dual bound.
This is a continuum inequality given exact weights and
integrals; the tabulated outer endpoints use numerical
quadrature and floating-point LP duals rather than certified
interval arithmetic. It bounds only the same frozen
kinematic model, not a time-dependent or cosmological result.

Use `early_pole_daughter_slope_limited_witness.py` for fixed
caps (options `--two-moments-only`, `--match-exterior` and
`--with-outer`), and
`early_pole_daughter_minimum_slope_threshold.py` at
`--subdivision 128` for the separate slope objectives.
Use the four matched frozen initial fixtures described in the
[prior envelope report](INDEPENDENT-EARLY-POLE-DAUGHTER-FROZEN-MOMENT-ENVELOPE-2026-11-27.md).

For hard depth four and 32 cold modes, the 512-subdivision
relaxation gives approximate paired upper endpoints
-0.004003 at (L=64), -0.001488 at (L=96), and +0.000360
at (L=128). Thus on this numerical integration the first two
exclude a positive initial rate, while the last is inconclusive.
The pointwise sign inclusion at much higher (L) follows from
actual admissible finite profiles, not from a positive outer
endpoint. The report retains the coarser relaxations, all tested
cap statuses and the weakest verified positive witness margin
(approximately (9.9993	imes10^{-6})).
