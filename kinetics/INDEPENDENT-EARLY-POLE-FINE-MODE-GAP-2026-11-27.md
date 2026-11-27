# The sampled late projection gap is coherent across cold modes

Replay the saved instantaneous source of each of the eight scalar momentum modes in the [matched fine five-channel 216-step pair](INDEPENDENT-EARLY-POLE-FINE-EXTENDED-MEMORY-2026-11-27.md). The archived-step branch exceeds the count-conserving continuous branch in all eight modes at each of the nine saved snapshots. Thus the small terminal total source gap is not due to cancellation across the sampled modes. It falls from .0217322 initially to 1.96071e−6 at step 8, then rebounds to .000182898 at step 216; the last gap and its absolute-bin sum coincide.

This differs from the [scalar self-redistribution test](INDEPENDENT-EARLY-POLE-SCREENED-SCALAR-MODE-CANCELLATION-2026-11-27.md), where large shifts in opposite directions canceled. The initial-state and scalar-self comparisons use different experiments and time horizons, so their magnitudes are not directly comparable. This checks nine stored snapshots at a single hard grid and step width, rather than every unsampled time, different projections, or the continuum.

From `kinetics/`, extract the two branch objects from `results/early-pole-joined-fine-216-memory-2026-11-27.json.gz` (keys `branches.steps` and `branches.cell-conservative`) into individual JSON files, then run:

```sh
OPENBLAS_NUM_THREADS=1 python early_pole_joined_mode_projection_gap_reference.py steps.json cell-conservative.json
```

Raw: `results/early-pole-joined-fine-216-mode-gap-2026-11-27.json`.
