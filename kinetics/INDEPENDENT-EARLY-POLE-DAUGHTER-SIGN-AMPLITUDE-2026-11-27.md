# Hard refinement raises the occupation change needed to reverse the initial sign

The [fixed-moment envelope](INDEPENDENT-EARLY-POLE-DAUGHTER-FROZEN-MOMENT-ENVELOPE-2026-11-27.md)
contains both signs when both daughter profiles can vary. To
measure how far such a witness is from the original fits,
minimize the maximum pointwise absolute change in occupation
inside the same physical interval. Keep each cell's own number,
bare energy and centered reciprocal-energy moment, with the same
frozen native packet geometry, parent and scalar seed.

In each small subcell impose the amplitude constraint using both
endpoints of the fitted occupation. The fitted logit is monotone
within a cell, so the endpoints bound the difference everywhere
between them for a piecewise-constant trial occupation. Require a
strictly positive total rate contrast, then check positivity,
moment residual and the actual amplitude after numerical moment
repair. These are finite witnesses, not exact continuum infima.

| Hard depth × cold modes | Both daughter profiles variable | Native only variable | Local only variable |
|---|---:|---:|---:|
| 1 × 8 | .01970 | .03926 | .04061 |
| 1 × 32 | .01970 | .03931 | .04055 |
| 4 × 8 | .13328 | Not feasible on this finite grid | Not feasible on this finite grid |
| 4 × 32 | .17615 | Not feasible on this finite grid | Not feasible on this finite grid |

Values are absolute occupations (for example, .020 is two
percentage points, not 2% relative to a background occupation).
The table uses 24 subcells per packet segment; the 12-subcell
witnesses give similar amplitudes. Independently, with the other
arm held at its original fit, the integrated numerical dual upper
values at hard depth four remain negative: native-only
-5.798e-3 / -6.882e-3 and local-only -5.621e-4 / -3.037e-3
at 8 / 32 cold modes. These outer estimates use ordinary
quadrature, not a certified error interval.

Thus substantially larger allowed occupation changes in both
arms are needed for the refined *frozen* sign witness. No result
shows that either altered shape occurs, that it is smooth,
or that it persists in time. Use
`early_pole_daughter_minimum_sign_amplitude.py` on the four
fixtures with `--subdivision 12` and 24. Combine the two outputs
and the earlier compressed envelope receipt with
`early_pole_daughter_minimum_sign_amplitude_compare.py`.
