# Hard-depth suppression is smaller in the second uninterrupted window

Continue the actual [depth-four cold32 pair](INDEPENDENT-EARLY-POLE-DAUGHTER-FINITE-HARD-COLD-CROSS-2026-11-27.md)
without resetting at the first window endpoint. Keep the same
native and local daughter initialization, both elastic channels,
32 cold modes, reporting spacing and tolerances. Compare to
the [depth-one double window](INDEPENDENT-EARLY-POLE-DAUGHTER-COLD32-DOUBLE-WINDOW-2026-11-27.md).
Depth four reproduces its independently published short prefix
exactly in the native arm and within 1.09e−19 per mode in the
local arm. The maximum double-window energy ledger defect
is 2.45e−15.

| Hard depth | First-window paired mean | Second-window paired mean | Full uninterrupted mean |
|---:|---:|---:|---:|
| 1 | +.00050819436 | +.00062641948 | +.00056730692 |
| 4 | +.00033440633 | +.00061307640 | +.00047374136 |

Depth-four suppression relative to depth one changes from 34.2%
in the first window to 2.13% in the second. The mixed finite
hard-depth and time correction is +.00016044494 in the same
rate units. All 32 depth-four net mode increments are positive
in each window and all sixteen second-window panel-by-panel
paired gains are positive. Nevertheless, each depth-four mode
contributes less paired gain than at depth one in both windows.
The cumulative sign change stays bracketed by [1.8750,2.0833]e−7
at depth four, later than [1.2500,1.4583]e−7 at depth one.

The depth-four-minus-depth-one contrast is negative on all but the
first panel of the complete run. In the second window its magnitude
decreases on every successive panel, from 2.61e−5 to 6.40e−6 in
the panel-mean rate. This is a finite-difference diagnostic, not an
instantaneous derivative.
 
The hard refinement acts differently in the two tested windows.
The small late correction does not establish hard-mesh convergence,
a universal late time, a physical expansion trajectory or a
prediction for the actual cosmological rate.

To regenerate the raw depth-four arms from the kinetics directory:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
i=results/early-pole-first-window-hard-depth4-cold32-initializer-2026-11-27.json.gz
a=(--duration 6.666666666666667e-7 --panels 32 --rtol 3e-9 --base-atol 2e-14 --record-mode-history)
python -W error early_pole_joined_simultaneous_chart_reference.py "$i" "${a[@]}"
python -W error early_pole_joined_simultaneous_chart_reference.py "$i" "${a[@]}" --local-from-paired-archive
```
For the recorded exact comparison:
```sh
OPENBLAS_NUM_THREADS=1 python -W error early_pole_daughter_two_window_compare.py \
 --hard-depth 4 \
 --short results/early-pole-unsplit-hard4-cold32-{native,local-daughter}-2026-11-27.json \
 --long results/early-pole-unsplit-hard4-cold32-double-window-{native,local-daughter}-2026-11-27.json.gz
OPENBLAS_NUM_THREADS=1 python -W error early_pole_daughter_hard_time_cross.py \
  results/early-pole-daughter-cold32-double-window-2026-11-27.json \
  results/early-pole-daughter-hard4-two-windows-2026-11-27.json
```
[Full crossed receipt](results/early-pole-daughter-hard-time-cross-2026-11-27.json).
