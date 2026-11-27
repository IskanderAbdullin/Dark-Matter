# The angularly refined paired gain persists after hard-depth refinement

Complete the missing corner of the
[depth-one angular evolution](INDEPENDENT-EARLY-POLE-DAUGHTER-ANGULAR78-HELD-CHART-2026-11-27.md).
On hard depth four, evolve both native and local hard-state controls
with the same hashed 78-mode angular mesh and exact cold seed used at
depth one. The local replacement spans eight nested cells over the
same daughter interval. The ordinary cold32 companions use the same
time mesh and solver settings.

| Hard depth | Paired mean, cold32 | Paired mean, angular78 | First cumulative sign-crossing bracket, both meshes |
|---:|---:|---:|---:|
| 1 | +.00050819435 | +.00050765147 | [1.250, 1.458]e−7 |
| 4 | +.00033440635 | +.00033400181 | [1.875, 2.083]e−7 |

The angular substitution reduces the hard4 mean by about 0.121%.
Hard-depth refinement reduces the angular mean by about 34.2%;
the extra mixed four-corner correction is +1.38337e−7.
All 78 hard4 terminal paired increments are positive, but every one
is smaller than at hard1. At this sampling neither angle substitution
moves the cumulative sign-crossing bracket. Across the eight archived
arms the largest energy and entropy ledger residuals are 2.89e−15 and
1.38e−14; the independent compensated scalar sums are within seven ulps.

These are actual paired first-window trajectories, not a bound on the
continuum hard or angular limit, a physical time of production, or a
cosmological abundance.

From the kinetics directory, generate the two missing trajectories
with `early_pole_daughter_angular_held_chart.py`, the published
`early-pole-first-window-hard-depth4-cold32-initializer-2026-11-27.json.gz`
and `early-pole-daughter-shared-cold-angular-sampling-2026-11-27.json`
under `results/`; add `--local` for the daughter arm.

Cross-check directly with the retained raw data:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error early_pole_daughter_angular_hard_cross.py \
  --hard1-angular results/early-pole-unsplit-hard1-angular78-{native,local-daughter}-2026-11-27.json.gz \
  --hard1-ordinary results/early-pole-unsplit-cold32-{native,local-daughter}-2026-11-27.json \
  --hard4-angular results/early-pole-unsplit-hard4-angular78-{native,local-daughter}-2026-11-27.json.gz \
  --hard4-ordinary results/early-pole-unsplit-hard4-cold32-{native,local-daughter}-2026-11-27.json \
  --initializer results/early-pole-daughter-angular-mode-initializer-2026-11-27.json
```
[Cross receipt](results/early-pole-daughter-angular-hard-cross-2026-11-27.json).
