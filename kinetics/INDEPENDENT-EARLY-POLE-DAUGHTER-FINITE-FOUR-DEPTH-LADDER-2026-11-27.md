# Actual depth eight on the matched 32-mode first window

Extend the [actual three-depth ladder](INDEPENDENT-EARLY-POLE-DAUGHTER-FINITE-THREE-DEPTH-LADDER-2026-11-27.md)
with two independently evolved depth-eight arms, native and local
daughter archive replacement. Hold the 32 cold momenta and initial
cold counts, 16 report panels, first exposure window, integration
tolerances and both elastic channels fixed. The replacement covers the
same physical daughter momentum interval, 1.4138888889 to 1.4666666667
in archive units, now divided into 16 daughter cells.

| Hard depth | Paired mean (local minus native) | Cumulative sign-crossing bracket |
|---:|---:|---:|
| 1 | +.00050819436 | [1.2500, 1.4583]e−7 |
| 2 | +.00045500547 | [1.4583, 1.6667]e−7 |
| 4 | +.00033440633 | [1.8750, 2.0833]e−7 |
| 8 | +.00031631524 | [1.8750, 2.0833]e−7 |

Depth 4→8 lowers the mean by 5.41%; this is about 0.150 times
the absolute depth 2→4 change. The earlier corrections grew
from 1→2 to 2→4, so this one smaller correction neither
establishes an asymptotic discretization law nor licenses
a continuum extrapolation. All 32 terminal mode effects still
decrease on every adjacent refinement and remain positive.
The sampled crossing bracket does not move at this last refinement.
Energy, entropy, and hard particle number defects are at most
2.0e−15, 7.8e−15, and 1.0e−15 across the two depth-eight arms.

This is only the **first window**. It cannot settle the late-window
small sign-changing depth effect seen in the
[three-depth joint time ladder](INDEPENDENT-EARLY-POLE-DAUGHTER-THREE-DEPTH-TIME-LADDER-2026-11-27.md).
The retained paired receipts are complete first-window trajectories,
not frozen static derivatives. The depth-eight time-extension control
must compare separately evolved arms over both windows.

With the retained initializer, each long replay starts from the same
32 cold modes:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
seed=results/early-pole-first-window-hard-depth8-cold32-initializer-2026-11-27.json.gz
python -W error early_pole_joined_simultaneous_chart_reference.py "$seed" --record-mode-history
python -W error early_pole_joined_simultaneous_chart_reference.py "$seed" --record-mode-history --local-from-paired-archive
```

For the quick four-depth check from retained raw:
```sh
python -W error early_pole_daughter_finite_three_depth_ladder.py \
  --hard1 results/early-pole-unsplit-cold32-{native,local-daughter}-2026-11-27.json \
  --hard2 results/early-pole-unsplit-hard2-cold32-{native,local-daughter}-2026-11-27.json.gz \
  --hard4 results/early-pole-unsplit-hard4-cold32-{native,local-daughter}-2026-11-27.json \
  --hard8 results/early-pole-unsplit-hard8-cold32-{native,local-daughter}-2026-11-27.json.gz
```
[Full matched result](results/early-pole-daughter-finite-four-depth-ladder-2026-11-27.json).
