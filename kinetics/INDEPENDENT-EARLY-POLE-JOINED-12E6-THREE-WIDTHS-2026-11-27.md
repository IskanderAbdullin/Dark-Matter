# The early count difference is concentrated before the rate histories nearly coincide

Match the first 1.2e−5 of the [216- and 432-step accepted long trajectories](INDEPENDENT-EARLY-POLE-JOINED-N432-TWO-HARD-INPUTS-2026-11-27.md) against new independent 144-step short runs at one-quarter the original width. Each lineage has the same initial hard moments, cold bins, and evaluated source per mode. Both short runs use the already verified 5e−15 elastic acceptance and separately pass terminal energy and entropy (energy drift below 2e−14). The two wider controls are prefixes of accepted longer trajectories with exactly the same initial data and collision controls, not forecasts.

| Width in exposure | Native integrated cold gain | Continuous integrated cold gain | Gap / native |
| --- | ---: | ---: | ---: |
| 3.3333e−7 | 1.08941067e−6 | 1.08553222e−6 | 0.3560% |
| 1.6667e−7 | 1.09016795e−6 | 1.08810837e−6 | 0.1889% |
| 8.3333e−8 | 1.09053443e−6 | 1.08938554e−6 | 0.1054% |

At the finest width the relative interval-rate gap stays below 0.1% after exposure 5.833e−7 within this window. The accumulated signed count gap is already 1.2205e−9 at the end of the first original coarse interval. It reaches a maximum 1.2776e−9 at exposure 3e−6 and falls back to 1.1489e−9 at 1.2e−5 as the continuous run gains at a slightly higher average rate. Nonetheless, the parent and daughter relative occupation L1 differences at the end remain around .66% and .30% over their full common hard supports. A nearly equal cold source does not establish equality of the spectra, freeze-out, or long-time convergence.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for p in steps cell-conservative; do
  python -W error early_pole_cold_screened_split_reference.py fine \
    --duration 1.2e-5 --count 144 --pair --same-soft-strength 1 --compton-strength 1 \
    --conversion-strength 1 --screening-update full --projection "$p" --snapshots \
    --elastic-energy-budget 5e-15 --checkpoint "$p-short-144.ckpt.json" > "$p-short-144.json"
done
python -W error early_pole_joined_three_widths_memory_reference.py \
  results/early-pole-joined-fine-216-memory-2026-11-27.json.gz \
  results/early-pole-joined-fine-432-time-refinement-2026-11-27.json.gz \
  steps-short-144.json cell-conservative-short-144.json three-widths.json.gz
```

Raw: `results/early-pole-joined-12e6-three-widths-2026-11-27.json.gz`.
