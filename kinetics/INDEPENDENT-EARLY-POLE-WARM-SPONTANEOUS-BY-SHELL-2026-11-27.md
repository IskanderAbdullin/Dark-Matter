# Most of the spontaneous warm quadrature bias comes from middle shells

Refine the [fixed-pole warm measure comparison](INDEPENDENT-EARLY-POLE-WARM-SPONTANEOUS-QUADRATURE-2026-11-27.md) in outer q³ quadrature, keeping radial order 8 and recording the eleven individual native shell contributions. Both hard reconstructions use their exact final n216 snapshot and its unchanged original scalar pole.

| Hard reconstruction | Marker-minus-resolved spontaneous rate, q12 | q24 | q48 |
| --- | ---: | ---: | ---: |
| Native steps | 1.3385038e−8 | 1.3390127e−8 | 1.3391521e−8 |
| Continuous | 1.3439612e−8 | 1.3444840e−8 | 1.3446272e−8 |

At q48 the five shells spanning q=.4348–2.7424 contribute 76.0–76.1% of the signed gap. The highest selected shell, with edges 3.9637–5.7288, contributes in the opposite direction, offsetting 5.6–5.7% of the net gap. The q24→q48 change is 5.3–5.5 ppm of the resolved warm rate, and the resolved spontaneous scalar-energy source is about 4.46% lower than its single-marker evaluation. These comparisons locate sensitivity to the integration measure, not a new channel or a physical temperature/abundance prediction.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for p in steps cell-conservative; do
  python -W error early_pole_warm_seed_measure_reference.py \
    results/early-pole-joined-fine-216-memory-2026-11-27.json.gz \
    --projection "$p" --outer-orders 12 24 48 --radial-orders 8
done
```

Raw: `results/early-pole-warm-spontaneous-by-shell-2026-11-27.json`.
