# Refine the first joined transient on its own time interval

The [full-window 216/432 comparison](INDEPENDENT-EARLY-POLE-JOINED-N432-TWO-HARD-INPUTS-2026-11-27.md) found that the first interval rate is less stable than the cumulative scalar gain. Start from exactly the same initial hard moments, cold mode counts, hard supports and active response. Restrict to (t=3.3333333333333335\times10^{-7}), the first original coarse width, and use 16 substeps with the stricter elastic energy budget. Compare interval averages, not the unequal instantaneous sources of the two inputs at the start.

| Steps per original interval | Archived-step hard input | Count-conserving continuous hard input |
| --- | ---: | ---: |
| 1 | .07378978827 | .06247869275 |
| 2 | .07010413607 | .06380618712 |
| 16 | .06694243427 | .06528878623 |

Refining from two to sixteen substeps moves the native average down 4.51% and the continuous average up 2.32%. The difference between them, measured relative to the native result, decreases from 15.33% to 8.98% to 2.47%. The two instantaneous archived sources at time zero (.07863619 and .05690395) were not changed. Both 16-step clipped histories pass the original terminal energy/entropy gates. This resolves how a coarse time interval can exaggerate the hard-input difference, but does not prove that the finest interval or an instantaneous source has converged. In particular it is not a full new evolutionary history.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for p in steps cell-conservative; do
  python -W error early_pole_cold_screened_split_reference.py fine \
    --duration 3.3333333333333335e-7 --count 16 --pair --same-soft-strength 1 --compton-strength 1 \
    --conversion-strength 1 --screening-update full --projection "$p" --snapshots \
    --elastic-energy-budget 5e-15 > "$p-first16.json"
  python -W error early_pole_joined_initial_window_refinement_reference.py \
    results/early-pole-joined-fine-216-memory-2026-11-27.json.gz \
    results/early-pole-joined-fine-432-time-refinement-2026-11-27.json.gz \
    "$p-first16.json" --projection "$p"
done
```

Raw: `results/early-pole-first-window-16-refinement-2026-11-27.json.gz`.
