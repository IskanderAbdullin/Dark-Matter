# How fragile is explicit-population CDF dominance to differing bin shapes?

Keep the original display counts from the three exact four-corner
rectangles. In every bin, let the same fraction of each corner's
measure have an independent conditional distribution, with the
remainder sharing a conditional distribution across all four corners.
Measure the CDF L1 amplitudes of the two signed effects independently
at unchanged counts.

The sum of absolute perturbations in the population and kernel CDFs
equals the larger of the diagonal and off-diagonal corner perturbation
differences. If two counts are a and b, the cube bound
`|a(x−z)−b(y−z)|<=max(a,b)` for conditional CDFs x,y,z gives a joint,
rather than separate, perturbation budget. Use this bound to certify a
sufficient release fraction for population dominance. For a converse
witness, construct a common endpoint mixture near the minimum population
CDF amplitude, then shift the explicitly released portions of the two
kernel groups towards opposite interior endpoints of each bin.

| Charts | Sufficient release bound (forward number) | Explicit reversal at |
|---|---:|---:|
| S2/S4 | 1.26597% | 1.30785% |
| S4/S8 | 0.336430% | 0.339210% |
| S8/S16 | 0.0849978% | 0.0852341% |

Thus a large dominance ratio in the full common-shape test is not an
unconditional attribution if the four in-bin shapes can differ. The
S8/S16 bracket uses the measured mesh-difference guard. It indicates
sensitivity under this mathematical family, not that numerical or
physical histories actually realize the adversarial shapes.

```sh
cd kinetics
d=2026-11-28
python -W error early_pole_endpoint_cold_spatial_shape_release.py --rectangles \
 "results/early-pole-s2-s4-population-edge-exact-rectangle-$d.json" \
 "results/early-pole-s4-s8-population-edge-exact-rectangle-$d.json" \
 "results/early-pole-cold-s8-s16-population-edge-exact-rectangle-$d.json" \
 --shapes "results/early-pole-cold-spatial-common-shape-geometry-$d.json"
```
