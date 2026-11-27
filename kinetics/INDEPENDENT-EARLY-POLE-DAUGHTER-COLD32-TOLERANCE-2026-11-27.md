# Sharper explicit integration on the same 32-mode cold mesh

Use the same original hard mesh, local replacement in the same
two daughter cells, fixed response and elastic channels, both
cold seeds, momentum nodes, time horizon and reported times.
Change only the DOP853 solver tolerances from relative 3e-9
and base absolute 2e-14 to relative 3e-10 and base absolute
2e-15. The charged logit and cold moment component scales retain
the original chart prescription. Check the original positivity,
number, energy and entropy conditions, then verify that individual
cold counts sum back to the saved scalar history within four
spacings of the scalar count.

Evaluate the sign terminally for each of the 32 output modes and
record the maximum tolerance sensitivity relative to the tightened
effect. Express the smallest effect in units of the corresponding
seed's stored floating-point spacing; this guards against a sign
reported only from a rounding artifact but does not bound spatial
or temporal discretization error by itself.

Use `early_pole_joined_simultaneous_chart_reference.py` on the
accepted original-hard 32-mode initial bundle with
`--record-mode-history`. Add `--local-from-paired-archive` for the
local arm, and `--rtol 3e-10 --base-atol 2e-15` for the tight pair.
Use `early_pole_daughter_cold32_tolerance_audit.py --base NATIVE LOCAL --tight TIGHT_NATIVE TIGHT_LOCAL`
for the held comparison.

## Result

All 32 terminal local-minus-native count shifts remain positive;
none changes any saved cumulative crossing bracket. The paired
mean changes from +5.08194357315e-4 to +5.08194354958e-4
(-4.64 parts per billion relative). The largest individual
terminal shift change relative to its tight magnitude is
2.61e-7. The smallest positive tight paired count is 7.3004e-21;
the least-resolved of the 32 shifts still exceeds 4.36e7 of
its initial mode's stored floating-point spacings. The counts
are therefore numerically resolved by this tolerance check;
this is not a proven bound for other tolerances or changed
hard support.

The native/local RHS calls increase from 1577/1541 to 2225/2177.
All four receipts reproduce the scalar sum within two stored
spacings. Both tighter raw files and the per-mode audit are
included in `results/`. Generate the original-hard 32-mode input
from the original eight-mode archive using the
[static rebin tool](early_pole_rebin_held_initial_cold_support.py).
