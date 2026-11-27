# Both elastic substep and outer split spacing affect the local daughter result

In the same held-response graph, set four backward elastic
substeps in each of the original 16 macro steps. Leave the slow
half steps unchanged. Compare native and two-cell daughter
preparations, then compare against full 64-step refinement with
one elastic step per macro interval. Both refined split variants
use 64 elastic solves with identical inner step width
(5.208333333333334e-9), but place the slow half steps differently.

| Macro × inner elastic steps | Mean cold-source shift, local minus native |
|---|---:|
| 16 × 1 | -1.00295605e-5 |
| 16 × 4 | +2.12757801e-4 |
| 64 × 1 | +3.84570917e-4 |
| Simultaneous held graph | +5.17245309e-4 |

Refining only the backward elastic step reverses the sign and
recovers 56.5% of the observed change from 16 × 1 to 64 × 1.
The remaining difference at the same inner step is 1.71813e-4
in the mean, with a maximum paired cumulative difference of
5.73e-11 on the saved grid. The 16 × 4 cumulative local effect
turns positive between 2.292e-7 and 2.500e-7; the 64 × 1 interval
is 1.667e-7 to 1.875e-7. Thus neither shortening only the inner
step nor treating outer interleaving as irrelevant is sufficient
for the tested benchmark. The percentages describe two numerical
controls, not physical shares or a two-variable error bound.

If a leading outer-width-squared term is valid at the fixed
inner width, the two-point cancellation
(16 Z(64 × 1) - Z(16 × 4)) / 15 gives +3.96025e-4.
Its difference from the simultaneous value is -1.21220e-4.
This is a conditional diagnostic, not an order measurement:
the hard remapping and finite slow solves can affect the outer
behavior. A 32 × 2 arm at the same inner width is a useful test
because the observed outer difference is substantial.

Both subcycled arms pass the terminal energy and entropy gates.
The largest recorded individual elastic defect is 4.89e-15 at a
5e-15 budget. An aggregate signed defect for the four substeps
within a macro step is stored separately, and is not an individual-step
budget measurement.

From `kinetics/`, reproduce the two new arms with
`early_pole_local_daughter_elastic_subcycle_reference.py` and
`--case native` or `--case local-daughter`. Compare against the
accepted matched split and simultaneous pairs with
`early_pole_daughter_subcycle_width_compare.py`.
