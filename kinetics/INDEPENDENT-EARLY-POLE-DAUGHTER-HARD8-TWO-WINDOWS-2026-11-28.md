# Two contiguous held windows at hard depth eight

Continue the archived [two-window comparison](INDEPENDENT-EARLY-POLE-DAUGHTER-HARD4-TWO-WINDOWS-2026-11-27.md) on the 32 cold-mode chart without restarting either distribution at the middle. Both full depth-eight integrations completed, with the same time tolerances and separate native and local-daughter occupations. The first window independently replays the published short depth-eight histories within 4.34e−19 per mode, under the same energy, particle-number and entropy gates. The two raw reports are retained.

| Hard depth | First-window paired mean | Second-window paired mean |
| --- | ---: | ---: |
| 1 | .0005081944 | .0006264195 |
| 2 | .0004550055 | .0006125300 |
| 4 | .0003344063 | .0006130764 |
| 8 | .0003163152 | .0006150038 |

The second-window mean changes by +1.9274e−6 from depth four to eight. The first mean changes by −1.8091e−5 and has a larger remaining depth dependence. At depth eight the paired total crosses zero in the same reported bracket as depth four, [1.875e−7,2.08333e−7]. These are finite held controls; they do not establish aperture convergence or supply a full cosmological time history.

Replay the [depth-eight paired receipt](results/early-pole-daughter-hard8-two-windows-2026-11-28.json) and the [four-depth receipt](results/early-pole-daughter-four-depth-time-ladder-2026-11-28.json) from `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
python -W error early_pole_daughter_four_depth_time_ladder.py \
  results/early-pole-daughter-cold32-double-window-2026-11-27.json \
  results/early-pole-daughter-hard2-two-windows-2026-11-27.json \
  results/early-pole-daughter-hard4-two-windows-2026-11-27.json \
  --short-8 results/early-pole-unsplit-hard8-cold32-{native,local-daughter}-2026-11-27.json.gz \
  --long-8 results/early-pole-unsplit-hard8-cold32-double-window-{native,local-daughter}-2026-11-28.json.gz
```
