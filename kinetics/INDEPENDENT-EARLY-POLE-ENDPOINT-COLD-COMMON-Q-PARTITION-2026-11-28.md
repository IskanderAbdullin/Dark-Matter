# Shared cold q boundaries, normalized profiles and conditional number neutrality

Use the complete, previously verified [S1/S2/S4 cold panels](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-S4-COLD-BAND-2026-11-28.md), preserving the fine charge-d2 chart and each already computed numerical pole. For each adjacent pair, use the intersection of their **exact** q boundaries; no profile is interpolated across a partial panel. Compare both accepted radial/angular mesh pairs using the same common boundaries. Each positive profile is normalized to its own band total before comparing shape.

| Comparison | Exact shared q intervals | Forward number TV lower bound | Inverse number TV lower bound | Forward excess TV lower bound |
| --- | ---: | ---: | ---: | ---: |
| S1→S2 | 10 | 8.502% | 8.579% | 9.617% |
| S2→S4 | 17 | 2.024% | 2.023% | 1.908% |

For the second pair alone, the q interval [3.30599941e−6,3.97454049e−6] loses 2.024 percentage points of the normalized gross forward share (from 47.1174% to 45.0934%), even though the raw band total changes −0.33018%. The different common partitions make these TV figures lower bounds for their own pairs, not proof of a convergence rate or a comparison of exact continuum TV. Similarly, quantile-bin support gives conservative ranges for W1 distance, not an intrapanel fit.

Next evaluate each chart at its own *conditional* single band-number-neutral occupation, (n_N=sum B/sum D). Its total normalized (n_N D-B) is zero, but common q groups have opposite signed net contributions. In the S2/S4 shared partition the positive-sum magnitudes are 0.38574% and 0.36975% of the corresponding gross inverse band integral; the negative sums precisely cancel them. This is a finite-panel lower estimate of the positive/negative redistribution at a diagnostic occupation. Neither "neutral" nor this cancellation is a measured photon population, detailed balance at every q, or a physical finite-band source. Individual fine/coarse mesh receipts and finite-group exact checks are provided.

Reproduction from `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_cold_q_partition
python -W error early_pole_endpoint_cold_q_partition.py \
 --left results/fine-charge2-scalar2-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --right results/fine-charge2-scalar4-endpoint-continuum-cold-{1,2}.jsonl.gz --levels 2 4
```
