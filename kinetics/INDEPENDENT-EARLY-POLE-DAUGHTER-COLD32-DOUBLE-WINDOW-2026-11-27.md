# The local cold gain persists in a second unreinitialized held interval

Continue both the native and two-cell daughter-substitution hard-depth-one
controls from the same prepared 32-bin cold states as the
[first held interval](INDEPENDENT-EARLY-POLE-UNSPLIT-LOCAL-COLD32-2026-11-27.md).
The new calculations run from exposure zero to 6.66667e−7, with 32 equally
spaced output panels and no reset at the original 3.33333e−7 endpoint.
The original 17 sampling times, exact seed and earlier per-mode
trajectories replay, with a maximum per-mode discrepancy below 4.34e−19
on either arm. The largest energy ledger residual is 1.12e−15.

| Exposure interval | Mean paired cold rate (local minus native) | Positive net paired modes |
|---|---:|---:|
| First interval | +.00050819436 | 32/32 |
| Second interval only | +.00062641948 | 32/32 |
| Whole uninterrupted run | +.00056730692 | 32/32 |

The second-window mean excess is 23.26% larger than the first-window
mean. The cumulative paired gain is positive at every reported sample
in the second interval and rises between every pair of its samples.
The only sampled negative-to-nonnegative crossing in the full trajectory
remains at [1.25000e−7,1.45833e−7]. These are statements about discrete
samples in the held model, not an instantaneous derivative bound,
a mapped FRW time, or indefinitely persistent production.

Replay from the kinetics directory. The initializer contains only the
exact prepared hard moments, common cold seed and supports of the two
original archives, not their later snapshots. The held solver checks
the native supports and the exact local replacement itself.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
init=results/early-pole-first-window-hard1-cold32-paired-initializer-2026-11-27.json.gz
solver=early_pole_joined_simultaneous_chart_reference.py
common=(--duration 6.666666666666667e-7 --panels 32 --rtol 3e-9 --base-atol 2e-14 --record-mode-history)
python -W error "$solver" "$init" "${common[@]}" > two-native.json
python -W error "$solver" "$init" "${common[@]}" --local-from-paired-archive > two-local.json
python -W error early_pole_daughter_two_window_compare.py \
  --short results/early-pole-unsplit-cold32-{native,local-daughter}-2026-11-27.json \
  --long two-native.json two-local.json
```

The [paired receipt](results/early-pole-daughter-cold32-double-window-2026-11-27.json)
also replays directly from the two compressed double-window trajectories.
