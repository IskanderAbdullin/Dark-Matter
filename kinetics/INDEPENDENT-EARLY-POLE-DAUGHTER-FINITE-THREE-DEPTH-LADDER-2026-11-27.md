# Actual depth two completes the hard-refinement ladder

Add the intermediate, uniformly subdivided daughter aperture to the
[actual depth-one and depth-four cold32 controls](INDEPENDENT-EARLY-POLE-DAUGHTER-FINITE-HARD-COLD-CROSS-2026-11-27.md).
All three levels use the identical 32-mode cold seed and momenta,
the same first time window, both elastic channels and matching
tolerances. The local replacement covers the identical physical
interval using two, four, or eight daughter cells, respectively.

| Hard depth | Paired mean (local minus native) | Cumulative sign-crossing bracket |
|---:|---:|---:|
| 1 | +.00050819436 | [1.2500, 1.4583]e−7 |
| 2 | +.00045500547 | [1.4583, 1.6667]e−7 |
| 4 | +.00033440633 | [1.8750, 2.0833]e−7 |

All 32 terminal paired mode gains remain positive on every mesh,
and all of them decrease at each adjacent hard refinement. The
magnitude of the second correction, depth 2→4, is about **2.27
times** that of 1→2. Thus the numerical sign is consistent across
these three levels, but the corrections do not shrink along the
ladder. These data do not justify a Richardson estimate or a
claim of an asymptotic hard-mesh limit. The observed crossing
also moves later with hard refinement; it is not a physical
time or a bound on the untested continuum.

To regenerate the hard2 preparations, start from the retained compact
paired eight-mode initializer and rebin only the cold coordinate:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
base=results/early-pole-first-window-hard-depth2-initializer-2026-11-27.json.gz
fine=results/early-pole-first-window-hard-depth2-cold32-initializer-2026-11-27.json.gz
python -W error early_pole_rebin_held_initial_cold_support.py "$base" "$fine"
python -W error early_pole_joined_simultaneous_chart_reference.py "$fine" --record-mode-history
python -W error early_pole_joined_simultaneous_chart_reference.py "$fine" --record-mode-history --local-from-paired-archive
```
For a direct comparison with the retained raw trajectories:
```sh
OPENBLAS_NUM_THREADS=1 python -W error early_pole_daughter_finite_three_depth_ladder.py \
  --hard1 results/early-pole-unsplit-cold32-{native,local-daughter}-2026-11-27.json \
  --hard2 results/early-pole-unsplit-hard2-cold32-{native,local-daughter}-2026-11-27.json.gz \
  --hard4 results/early-pole-unsplit-hard4-cold32-{native,local-daughter}-2026-11-27.json
```
[Full three-depth result](results/early-pole-daughter-finite-three-depth-ladder-2026-11-27.json).
