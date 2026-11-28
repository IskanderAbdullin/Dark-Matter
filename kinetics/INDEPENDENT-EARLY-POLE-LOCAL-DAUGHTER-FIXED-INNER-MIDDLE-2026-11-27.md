# Three outer widths at one identical elastic width

Use the [same held graph and two-cell prepared control](INDEPENDENT-EARLY-POLE-LOCAL-DAUGHTER-ELASTIC-SUBCYCLE-2026-11-27.md).
Insert the intermediate 32 macro by 2 elastic profile between
16 by 4 and 64 by 1. This changes only where the slow half steps
are interleaved with the 64 backward elastic solves. The
comparison validates every individual elastic width from the
ledgers, the exact starting hard/cold data and the exposure grid.

| Outer macro steps × inner substeps | Local minus native mean |
|---|---:|
| 16 × 4 | +2.12757800960e-4 |
| 32 × 2 | +3.46611354202e-4 |
| 64 × 1 | +3.84570916778e-4 |
| Simultaneous, no splitting | +5.17245309084e-4 |

The ratio (16 × 4 minus 64 × 1)/(32 × 2 minus 64 × 1)
is 4.526; a single outer-width-squared term would predict 5.
The two outer-quadratic cancellations (using the 16-to-64 and
32-to-64 pairs, respectively) give +3.96025e-4 and +3.97224e-4,
different by 1.20e-6. Both remain about 1.20e-4 below the
simultaneous result. Thus interleaving matters, and changing it
alone at this inner width does not remove the remaining error.
The ratio is a finite diagnostic rather than demonstrated order;
mapping thresholds and approximate slow solves can affect it.

The cumulative sign crossing of the intermediate control falls
between exposures 1.875e-7 and 2.083e-7, between the two other
split crossings. Largest individual intermediate elastic energy
defect is 4.663e-15 under the 5e-15 budget; original terminal
gates are also met. The specific next control changes only the
inner elastic width at 64 outer steps.

Reproduce with the previous subcycle command and
`--macro-steps 32 --elastic-substeps 2`. The validator's optional
`--middle` pair takes the native and local intermediate raw
outputs.
