# Continuous occupations inside the aperture still permit both frozen signs

Repeat the [frozen moment envelope](INDEPENDENT-EARLY-POLE-DAUGHTER-FROZEN-MOMENT-ENVELOPE-2026-11-27.md)
using continuous piecewise-linear daughter profiles instead of
discontinuous steps. Share the same nodal value across adjacent
subcells and across each neighboring hard-cell boundary within
the selected interval. Enforce each original count, energy and
centered reciprocal moment as before. Positivity holds between
nodes by affine interpolation and is checked again after correcting
the numerical moment residual.

At 24 subdivisions per allowed packet segment, the feasible
paired ranges are [-.088438,+.060998] and [-.089857,+.058752]
for hard depth one (8 and 32 cold modes), and
[-.027474,+.004981] and [-.025830,+.002233] at depth four.
Both signs also occur at 12 subdivisions, and the largest
normalized moment defect in the checked range is below 2.78e-17.

The witnesses are continuous within the selected interval.
They need not match outside it, have a bounded physical slope,
be differentiable at nodes, or arise dynamically. This removes
one numerical discontinuity objection only. Run
`early_pole_daughter_continuous_sign_witness.py` on the same
four fixtures at `--subdivision 12` and 24.
