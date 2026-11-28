# How large is the explicit-population shape effect with unknown shared bin shapes?

Use the exact four-corner physical-input substitutions on the same
original 17-bin display for S2/S4, S4/S8 and S8/S16. Within each bin
leave the conditional distribution unspecified, but require that it
be shared by the four counterfactual profiles. This deliberately
isolates differences in the measured bin counts; it does not cover
independent unknown shapes for the four corners.

For a signed bin-count difference, let the differences in the cumulative
sums at the two bin edges be a and b. The CDF inside is a+(b−a)F,
where F is the common conditional CDF. The smallest possible magnitude
is zero if the signs differ, otherwise min(|a|,|b|); the largest is
max(|a|,|b|). Multiply by bin width and sum. Concentration at edges, or
its interior limiting version, attains the bounds separately for each
effect. The total native contrast has the usual W1 interpretation;
the symmetric population and combined-kernel terms are signed CDF
contributions. Separate extrema need not share the same witness.

| Adjacent charts | Native forward W1 bracket | Lower ratio of population to combined-kernel CDF amplitude |
|---|---:|---:|
| S2/S4 | [8.7006e−9, 3.0013e−8] | 16,890 |
| S4/S8 | [2.2776e−9, 7.0943e−9] | 6,862 |
| S8/S16 | [5.7345e−10, 1.7514e−9] | 27,738 |

The ratios use the opposite extremal sides and conservatively account
for the measured two-mesh profile difference and a floating guard.
They compare resolved shape amplitudes under a common-bin convention
only. They do not identify actual intrabin transport or a physical
convergence rate.

```sh
cd kinetics
d=2026-11-28
python -W error early_pole_endpoint_cold_spatial_common_shape_geometry.py --rectangles \
 "results/early-pole-s2-s4-population-edge-exact-rectangle-$d.json" \
 "results/early-pole-s4-s8-population-edge-exact-rectangle-$d.json" \
 "results/early-pole-cold-s8-s16-population-edge-exact-rectangle-$d.json"
```
