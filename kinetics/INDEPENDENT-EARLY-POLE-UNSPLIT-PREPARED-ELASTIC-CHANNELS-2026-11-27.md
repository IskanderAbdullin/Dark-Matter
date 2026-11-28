# Elastic channel selection for the local and outside compensation responses

In the original held-response simultaneous chart, select both
elastic terms, just unlike scattering, just same-flavor scattering,
or neither. In every profile hold the native rate calculation, the
local aperture swap and the compensation-only input fixed apart
from that switch; keep the cold seed, geometry and reporting times
identical. The preflight checks that at each of the three initial
states summing the two instantaneous elastic vector fields matches
the full elastic field to less than 1e-16 relative, with no change
to the nonelastic right-hand side when unlike strength is toggled.

| Elastic profile | Local minus native mean | Complement only minus native mean |
|---|---:|---:|
| Neither | -1.55091162e-2 | -4.20642687e-7 |
| Unlike only | +5.11567184e-4 | -4.36550941e-4 |
| Same-flavor only | -5.71741699e-3 | -1.15018162e-4 |
| Both | +5.17245309e-4 | -4.28892890e-4 |

Unlike-only retains 98.9% of the joined positive local shift
and 101.8% of the joined complement suppression in this switch
test. Same-flavor-only still gives a negative local shift.
The finite both-minus-unlike-minus-same-plus-neither contrast is
-9.78602e-3 for local and +1.22256e-4 for the complement.
Consequently treating each channel's finite change relative
to the nonelastic profile as an additive contribution would be
incorrect. All twelve cases pass the original gates, with
maximum absolute energy defect about 1.6e-15.

Additivity of initial RHS vectors does not imply additivity of
evolved preparation contrasts. The comparisons are controlled
model switches, not unique assignments of physical production
to scattering channels or uncertainties for a physical relic.

Reproduce using the original fine paired archive and
`early_pole_joined_simultaneous_chart_reference.py` with either
`--elastic-channel unlike`, `--elastic-channel same`, both by
default, or `--without-elastic`. Supply accepted native, local and
compensation-only raw paths in that order for each channel set to
`early_pole_unsplit_prepared_elastic_channels_compare.py`.
