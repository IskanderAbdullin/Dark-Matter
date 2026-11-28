# Longer fine history with twice the time resolution

Repeat the [2.4e−5 five-channel fixed-scale comparison](INDEPENDENT-EARLY-POLE-JOINED-LOCAL-SPECTRAL-MEMORY-2026-11-27.md) with 144 rather than 72 macro steps. Both branches use the same fine hard momentum grid. The scalar initial counts, collision supports, hard channels and like-flavor frozen response are unchanged, and the unlike soft response is still refreshed. The comparison varies only the macro time grid and the original two initial projection choices.

| Initial reconstruction | Gain at 72 steps | Gain at 144 steps | Relative change |
| --- | ---: | ---: | ---: |
| Archived steps | 2.35740053e−6 | 2.35889018e−6 | +0.0632% |
| Count-conserving continuous | 2.35307536e−6 | 2.35638304e−6 | +0.1406% |

The difference between the integrated gains shrinks from 0.1835% to 0.1063% of the step result. The final interval rate difference stays 0.070%, and the comparable instantaneous final rates are .10884347 and .10876715. After exposure 6.67e−7 the 144-step interval averages all differ by less than 0.1% of the step rate. All reported interval gains remain positive.

The final relative occupation L1 remains 0.0871% and 0.1637% for the parent and daughter within [1.39,1.48]. On the full common live support, the differences remain 0.680% and 0.411%. Thus finer time steps refine the accumulated gain but do not turn a local spectral agreement into global equality. The energy residual is below 1.6e−13, and entropy mismatch below 1.6e−12. This is convergence evidence at fixed hard mesh, finite channels and short late exposure, not evidence for continuum convergence, loss of all initial information, or freeze-out.

From `kinetics/`, with one BLAS thread:

```sh
for p in steps cell-conservative; do
  OPENBLAS_NUM_THREADS=1 python early_pole_cold_screened_split_reference.py fine \
    --duration 2.4e-5 --count 144 --same-soft-strength 1 --compton-strength 1 \
    --conversion-strength 1 --screening-update full --pair --projection "$p" --snapshots > "$p.json"
done
OPENBLAS_NUM_THREADS=1 python early_pole_joined_history_memory_reference.py steps.json cell-conservative.json
OPENBLAS_NUM_THREADS=1 python early_pole_joined_snapshot_aperture_reference.py steps.json cell-conservative.json
```

Raw: `results/early-pole-joined-fine-144-memory-2026-11-27.json.gz`.
