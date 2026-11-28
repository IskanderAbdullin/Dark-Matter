# Live warm-marker redshift remains a small matched correction at twelvefold exposure

Continue the [short prescribed-scale held-polarization experiment](INDEPENDENT-EARLY-POLE-LIVE-WARM-REDSHIFT-2026-11-27.md) from the same fine final n216 snapshot. Keep the six-step short run's macro width, using 72 steps to exposure 2.4e−5. Match fixed vs prescribed log(a)=t with the same selected warm markers, radiation-only collision clock, held dielectric and other joined hard channels. This increases the exposure, not the resolution at the longer endpoint.

| Clock path | Warm gain | Cold gain (live warm arm) |
| --- | ---: | ---: |
| Fixed scale | 6.59757932979e−12 | 2.84561293605e−6 |
| Prescribed redshift | 6.59781697751e−12 | 2.84571514334e−6 |

The expansion path increases the warm gain by 36.0204 ppm and the cold gain by 35.9175 ppm relative to their respective fixed-scale results. The warm increment is close to twelve times the 3.0017-ppm short comparison. In the expanding arm the independently integrated hard redshift trace is 9.365e−9; subtracting it leaves an energy defect near 1.825e−12, below the unchanged 3e−12 gate. The fixed arm has a similar defect. The accepted short-time halving comparison does not establish time convergence at the longer endpoint. Neither path solves for the Friedmann clock nor updates the dielectric tables, and this comparison does not infer a late abundance or measure a tiny warm-to-cold effect from two separately solved populations.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
p=early_pole_live_scalar_exterior_frw_reference.py
python -W error "$p" "$a" --duration 2.4e-5 --panels 72
python -W error "$p" "$a" --duration 2.4e-5 --panels 72 --expands
```

Raw: `results/early-pole-live-warm-redshift-longer-2026-11-27.json`.
