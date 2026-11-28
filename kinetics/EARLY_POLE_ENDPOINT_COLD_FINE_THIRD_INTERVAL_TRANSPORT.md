# Native fine/third interval transport measures displacement separately from shape difference

Use the complete fine and third D2 cold panel ledgers and keep their respective normalized masses unchanged. An interval coupling charges only the distance between the two corresponding q supports, allowing arbitrary locations inside each panel. This finite problem is independently bracketed with a repaired feasible coupling and constrained potentials. An independent continuous CDF-envelope lower value confirms the same scale without solving the coupling problem.

| Cold gross response | Minimal compatible normalized TV | Support-only interval cost in q | Integrated CDF-envelope lower value |
| --- | ---: | ---: | ---: |
| Forward number | 26.23343% | 6.382744e−8 | 6.382744e−8 |
| Inverse number | 26.40743% | 6.421399e−8 | 6.421399e−8 |
| Forward excess | 31.49369% | 7.878841e−8 | 7.878841e−8 |

The forward interval cost is about 0.456% of the full width of the retained cold band, even though minimizing the mass difference itself still leaves a 26.2% TV gap. The two metrics measure different things. Neither establishes the true location inside the panels, treats the coarse and fine chart as the same underlying physical state, or measures a cosmological transfer.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results; d=2026-11-28
python -W error early_pole_endpoint_cold_fine_third_interval_transport.py \
 --fine "$p/fine-d2-endpoint-continuum-cold"-{1,2}.jsonl.gz \
 --third "$p/third-d2-endpoint-continuum-cold"-{1,2}.jsonl.gz \
 --states "$p/early-pole-cell-endpoint"-{fine,third}-d2-x05-state-$d.json.gz \
 --archive "$p/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json" \
 --fixture "$p/early-pole-cold-fine-third-d2-native-cdf-$d.json"
```
