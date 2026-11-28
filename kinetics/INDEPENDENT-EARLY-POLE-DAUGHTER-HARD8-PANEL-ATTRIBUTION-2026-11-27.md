# The depth-eight correction reverses within the first window

The [matched depth-four to depth-eight result](INDEPENDENT-EARLY-POLE-DAUGHTER-FINITE-FOUR-DEPTH-LADDER-2026-11-27.md)
lowers the first-window total. Resolve when that difference is
accumulated, reconstructing each paired difference from the four
raw mode histories on the same 32 modes and 16 report panels.

The first ten panel increments in the H8−H4 contrast are negative,
while the final six are positive. The sign changes between panel
midpoints 1.97917e−7 and 2.18750e−7. At the first positive total
panel increment, 17 of 32 mode increments have reversed. All 32
increase in each of the final five panels. This concerns the
*time increment* in the difference between depths, not the sign
of the already accumulated paired response.

The parts of the H8−H4 count accumulated in each equal quarter are
−5.32237e−12, −8.28230e−13, −2.48156e−14 and +1.45052e−13.
The first quarter accounts for about 88.3% of the final net reduction.
The cumulative H8 result therefore still lies below H4 at the end.
The terminal per-mode paired H8/H4 ratios range from .81245 to .97950.

This locates the phase of the numerical refinement effect in the
first window. A later positive increment does not establish the
second-window average or a continuum result. Both second-window
arms must be solved over that longer exposure with the same seeds.

Replay from kinetics with the retained raw data:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error early_pole_daughter_hard8_panel_attribution.py \
  --hard4 results/early-pole-unsplit-hard4-cold32-{native,local-daughter}-2026-11-27.json \
  --hard8 results/early-pole-unsplit-hard8-cold32-{native,local-daughter}-2026-11-27.json.gz \
  --ladder results/early-pole-daughter-finite-four-depth-ladder-2026-11-27.json
```
[Raw modal panel comparison](results/early-pole-daughter-hard8-panel-attribution-2026-11-27.json).
