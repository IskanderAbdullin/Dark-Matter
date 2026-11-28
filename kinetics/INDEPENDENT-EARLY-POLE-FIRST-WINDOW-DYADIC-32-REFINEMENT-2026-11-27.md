# First-interval mean approaches a much smaller hard-input gap

Refine the [first coarse interval](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-16-REFINEMENT-2026-11-27.md), duration 3.3333333333333335e−7, using the same prepared supports and full unlike response. Compare conservative accepted histories at 8, 16 and 32 substeps, with absolute elastic energy acceptance 5e−15. The starting hard moments, scalar bins and evaluated source per mode are exactly the same within each input lineage.

| Substeps on common interval | Native archived steps | Count-conserving continuous | Gap / native |
| --- | ---: | ---: | ---: |
| 1 | .07378978827 | .06247869275 | 15.33% |
| 2 | .07010413607 | .06380618712 | 8.98% |
| 8 | .06738090709 | .06505519482 | 3.45% |
| 16 | .06694243427 | .06528878623 | 2.47% |
| 32 | .06672255625 | .06540677346 | 1.97% |

The signed 8→16 and 16→32 differences nearly halve in both branches. The observed powers from those three values are .996 and .985. A first-order extrapolation gives .06650268 and .06552476 and a gap of about 1.47%. This extrapolation is a diagnostic, not a validated exact limit. The local starting instantaneous rates are not estimated by this calculation, and the mean across the first interval is not a late-time rate or a cosmological abundance.

An attempted native 16-step midpoint alternative did not pass its nonlinear acceptance gate on the first macro step (residual approximately 4.38e−11 versus the unchanged 5e−13 threshold). It is not used in this table and does not indicate a physical breakdown. The table concerns the tested backward-elastic sequence.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for p in steps cell-conservative; do
  for n in 8 16 32; do
    python -W error early_pole_cold_screened_split_reference.py fine \
      --duration 3.3333333333333335e-7 --count "$n" --pair --same-soft-strength 1 --compton-strength 1 \
      --conversion-strength 1 --screening-update full --projection "$p" --snapshots \
      --elastic-energy-budget 5e-15 > "$p-first-$n.json"
  done
  python -W error early_pole_first_window_cauchy_reference.py \
    results/early-pole-joined-fine-216-memory-2026-11-27.json.gz \
    results/early-pole-joined-fine-432-time-refinement-2026-11-27.json.gz \
    "$p-first-8.json" "$p-first-16.json" "$p-first-32.json" --projection "$p"
done
```

Raw: `results/early-pole-first-window-dyadic-32-refinement-2026-11-27.json.gz`.
