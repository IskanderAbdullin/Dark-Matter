# Separate the individual hard-refinement changes from their paired residual

The [first-window panel comparison](INDEPENDENT-EARLY-POLE-DAUGHTER-HARD8-PANEL-ATTRIBUTION-2026-11-27.md)
found the H8−H4 change in the local-minus-native rate flips sign
late in the window. Reconstruct both individual H8−H4 histories
from the four raw mode trajectories, before taking their difference.
All four use the same 32 scalar seeds, report times, controls and
conservation gates.

The terminal count changes are −7.19325e−11 for the native arm
and −7.79629e−11 for the local daughter arm; all 32 terminal modes
decrease in each arm. Taking their difference cancels 95.98% of
the sum of the absolute arm changes and leaves −6.03036e−12.
At the first panel the two signed arm mean corrections are
positive, but the local correction is smaller, so their paired
difference is already negative. After that first panel both
individual corrections are negative in every saved panel.

In the final panel the individual rate corrections are
−2.62825e−4 and −2.60814e−4, respectively. Their difference is
+2.01074e−6. The common negative rate change is about 130 times
as large as the paired differential change. The late paired
reversal therefore does not mean either individual arm reversed.

The tool reads the raw exposure length, so the same accounting
can be applied to the two separately evolved longer windows when
both are available. No particular second-window sign or convergence
order follows from this first-window result.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error early_pole_individual_arm_refinement_reference.py \
  --low results/early-pole-unsplit-hard4-cold32-{native,local-daughter}-2026-11-27.json \
  --high results/early-pole-unsplit-hard8-cold32-{native,local-daughter}-2026-11-27.json.gz \
  --paired-receipt results/early-pole-daughter-hard8-panel-attribution-2026-11-27.json
```
[Raw separate-arm receipt](results/early-pole-daughter-h8-individual-armed-refinement-2026-11-27.json).
