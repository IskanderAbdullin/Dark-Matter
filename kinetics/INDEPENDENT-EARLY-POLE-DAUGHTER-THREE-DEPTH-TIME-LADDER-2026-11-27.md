# The hard-depth correction reverses between the two held windows

Add a matched H2 native/local daughter pair over the same extended
exposure as the [existing H1/H4 two-window comparison](INDEPENDENT-EARLY-POLE-DAUGHTER-HARD-TIME-CROSS-2026-11-27.md).
The independently published [H2 short run](INDEPENDENT-EARLY-POLE-DAUGHTER-FINITE-THREE-DEPTH-LADDER-2026-11-27.md)
is recovered within 4.34e−19 per mode in each arm. The comparison
retains 32 cold modes, both elastic channels, the original time
slices, the initial state and the same tolerances. Each long run
continues without a reset after the first window.

| Hard depth | First-window paired mean | Second-window paired mean |
|---:|---:|---:|
| 1 | 0.000508194356 | 0.000626419484 |
| 2 | 0.000455005471 | 0.000612530033 |
| 4 | 0.000334406331 | 0.000613076397 |

At the first window both refinements reduce the paired gain; the
H2→H4 reduction is 2.27 times the H1→H2 reduction. In the second
window the H2→H4 correction turns positive (+5.46e−7), is only
3.93% of the magnitude of H1→H2, and increases all 32 mode
contributions (0.086–0.107% per mode). Thirteen of sixteen
second-window panel increments also increase from H2 to H4.
All represented modes still contribute positive net paired gains at
all three depths in each window. H2's cumulative sign crossing
falls in [1.4583,1.6667]e−7, between the reported H1/H4 intervals.
The H2 energy defects are at most 1.34e−15.

The later correction disagrees in sign with monotonic first-window
suppression: comparing only H1 to H4 hides the rebound. Neither
the small late difference nor the sign of the early sequence
establishes a continuum rate, converged mesh order, or a physical
expanding solution.

To regenerate each H2 long arm from kinetics:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
i=results/early-pole-first-window-hard-depth2-cold32-initializer-2026-11-27.json.gz
a=(--duration 6.666666666666667e-7 --panels 32 --rtol 3e-9 --base-atol 2e-14 --record-mode-history)
python -W error early_pole_joined_simultaneous_chart_reference.py "$i" "${a[@]}"
python -W error early_pole_joined_simultaneous_chart_reference.py "$i" "${a[@]}" --local-from-paired-archive
```
The recorded comparison is directly replayable:
```sh
OPENBLAS_NUM_THREADS=1 python -W error early_pole_daughter_two_window_compare.py \
  --hard-depth 2 \
  --short results/early-pole-unsplit-hard2-cold32-{native,local-daughter}-2026-11-27.json.gz \
  --long results/early-pole-unsplit-hard2-cold32-double-window-{native,local-daughter}-2026-11-27.json.gz
OPENBLAS_NUM_THREADS=1 python -W error early_pole_daughter_three_depth_time_ladder.py \
  results/early-pole-daughter-{cold32-double-window,hard2-two-windows,hard4-two-windows}-2026-11-27.json
```
[Raw three-depth receipt](results/early-pole-daughter-three-depth-time-ladder-2026-11-27.json).
