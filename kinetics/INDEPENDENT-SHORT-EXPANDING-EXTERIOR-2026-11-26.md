# Short expanding embedded-exterior control

Continue the original q16 and embedded q32 vacuum-source archives from
the common 12 MeV endpoint to the marked scale factor b=1.1066819197.
The q32 case preserves the original inner edges and adds four populated
exterior shells. Both cases retain the nonlinear below-pole elastic
operator and the retained on-shell parent, include radiation plus
retained dark gravity, and switch off late visible sources. The q16
and q32 archives are separately generated; this is not a paired
experiment that changes only the *later* boundary condition.

Use the compact cache/elastic implementation as an optional numerical
path. This does not change the backward-Euler residual criterion,
cutoffs, or relaxation law. It is possible to switch off elastic
scattering or use the independent parent-local pivot, but the table
below uses the compact elastic path and the previous prompt center.
The old q16 edge is an internal diagnostic on the wide run.

| Late intervals | q16 T_phi (MeV) | q32 T_phi (MeV) | q32 scalar number outside old edge |
| ---: | ---: | ---: | ---: |
| 1 | 9.97828 | 9.97921 | 1.91053e-5 |
| 2 | 9.94499 | 9.94546 | 1.81751e-5 |

The q32 scalar exterior starts at 2.21350e-5. The scalar temperature
difference from widening the grid is smaller than the time-step
dependence in this short interval. Neither exterior depletion nor
temperature at this short endpoint estimates the final relic density.
The total expanding count and pressure-energy ledgers are checked
separately from the exterior diagnostic.

## Reproduce

Prepare [both EOS inputs](INDEPENDENT-BENCHMARK-SETUP.md), then use
`12-vacuum-24` or `wide32-embedded-vacuum` as the archive stem in
`kinetics/results`:

```sh
export OPENBLAS_NUM_THREADS=1
python -W error late_archived_shape_continuation_reference.py EOS.dat results/"$STEM"-2026-11-26.json results/"$STEM"-shape-2026-11-26.json --steps 2 --stop 1.1066819197 --prompt-until 20 --radiation-table EOS.late.dat --dark-gravity --compact-elastic
```

The [fixed-time prompt transient](INDEPENDENT-PROMPT-EDGE-DYNAMIC-TRANSIENT-2026-11-26.md)
also passes a tighter Radau tolerance at the same endpoint, with the
exterior scalar fraction changing by less than 5e-18 in absolute
number-fraction units. This is independent of the expansion-step comparison.
