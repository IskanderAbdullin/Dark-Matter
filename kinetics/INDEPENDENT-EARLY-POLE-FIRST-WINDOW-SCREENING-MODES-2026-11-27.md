# Updating screening is a sub-ppm correction on the first joined interval

On the [matched first 16-substep interval](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-16-REFINEMENT-2026-11-27.md) keep both elastic channels present. In the unlike channel compare held response, an update of the static Debye normalization only, and the full stored response update. Hold the same-flavor tables in all three cases. The prepared state, scalar bins, hard geometry, strengths, time grid and 5e−15 elastic step acceptance are otherwise unchanged within each lineage.

| Hard input | Held mean cold rate | Static-only gain shift | Full gain shift |
| --- | ---: | ---: | ---: |
| Native steps | .0669424706293 | −0.542943 ppm | −0.543098 ppm |
| Count-conserving continuous | .0652888243885 | −0.584434 ppm | −0.584434 ppm |

The full-minus-static correction is −1.55e−10 of the held gain in the native input and unresolved in the count arithmetic for the continuous input. The maximum absolute fractional Debye-square change is about 3.9e−4; the maximum absolute shift among the reported soft packet-weight quantiles is about 8.6e−4. Those response shifts need not equal the integrated cold-source correction. All three controls pass the original terminal energy and entropy checks.

This finite early comparison separates the [large effect of switching on elastic collisions](INDEPENDENT-EARLY-POLE-FIRST16-ELASTIC-FACTORIAL-2026-11-27.md) from the much smaller effect of updating the tables while collisions remain on. No inference follows about retarded screening later in the trajectory, a full evolving plasma, or the final cosmological abundance.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for p in steps cell-conservative; do
  for s in frozen debye full; do
    python -W error early_pole_cold_screened_split_reference.py fine \
      --duration 3.3333333333333335e-7 --count 16 --pair --same-soft-strength 1 --compton-strength 1 \
      --conversion-strength 1 --projection "$p" --screening-update "$s" \
      --snapshots --elastic-energy-budget 5e-15 > "$p-$s.json"
  done
  python -W error early_pole_initial_screening_mode_reference.py --projection "$p" \
    --frozen "$p-frozen.json" --debye "$p-debye.json" --full "$p-full.json"
done
```

Raw: `results/early-pole-initial-screening-mode-2026-11-27.json`.
