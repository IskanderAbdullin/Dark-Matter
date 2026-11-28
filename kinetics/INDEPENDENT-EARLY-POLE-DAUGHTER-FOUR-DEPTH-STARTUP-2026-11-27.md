# Initialized rates do not predict the finite matched depth-eight response

On the [actual cold32 first-window ladder](INDEPENDENT-EARLY-POLE-DAUGHTER-FINITE-FOUR-DEPTH-LADDER-2026-11-27.md),
reconstruct the initial native and aperture-local daughter states from
the retained paired initializers. Independently evaluate the cold
right-hand side without taking a time step. Retain the identical 32
cold bins, the exact archive-matched hard cells, and only the intended
local replacement on the same physical interval. The depth-eight
initialized local rate exactly matches the independently stored
finite trajectory's initial-rate receipt.

| Hard depth | Initial local-minus-native rate | Actual first-window mean |
|---:|---:|---:|
| 1 | −.01625060 | +.00050819436 |
| 2 | −.01197832 | +.00045500547 |
| 4 | −.01223583 | +.00033440633 |
| 8 | −.01221517 | +.00031631524 |

All 32 initialized modal contrasts are negative at each depth,
while every terminal paired effect is positive. In the last
refinement the initial contrast becomes slightly less negative,
by approximately 2.07e−5, while the finite-window mean decreases
by 1.81e−5. The static startup trend cannot be used as a proxy for
even the sign of that temporal refinement.

At depth eight the highest eight momentum nodes contribute about
1.60% of the positive terminal count difference, versus 1.57% at
depth four. The paired-count weighted root-mean-square momentum
rises from 3.999e−6 to 4.014e−6 in the archived T units; it is a
kinematic moment of the positive differences, not the full evolved
distribution or a cosmological transfer prediction.

Replay from kinetics with the retained initializers and finite receipts:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error early_pole_daughter_four_depth_startup_control.py \
  results/early-pole-first-window-hard1-cold32-paired-initializer-2026-11-27.json.gz \
  results/early-pole-first-window-hard-depth{2,4,8}-cold32-initializer-2026-11-27.json.gz \
  --finite results/early-pole-daughter-finite-four-depth-ladder-2026-11-27.json \
  --local8 results/early-pole-unsplit-hard8-cold32-local-daughter-2026-11-27.json.gz
```
[Raw initialized and finite comparison](results/early-pole-daughter-four-depth-startup-control-2026-11-27.json).
