# Refined warm-marker redshift with explicit elastic energy closure

Repeat the [long prescribed expansion test](INDEPENDENT-EARLY-POLE-LIVE-WARM-REDSHIFT-LONGER-2026-11-27.md)
in the same fixed and expanding geometries, from the same archived fine final
hard snapshot and through the same exposure of 2.4e-5. The revised numerical
profile uses 144 rather than 72 macro steps and imposes a 5e-15 local
elastic energy gate. Both changes are applied together, so this is not a
separate estimate of either time discretization or solver energy error.

The independent large redshift term is separated from terminal conservation.
The four strict redshift-subtracted defects are all below 4.5e-15,
compared with roughly 1.825e-12 in the original runs. The maximum
recorded local elastic energy defect is 1.56e-15; no substep splits
were needed. The absolute hard count defects fall from roughly 4.55e-13
to below 8e-16.

| Numerical profile | Expanding minus fixed, warm gain | Expanding minus fixed, cold gain |
|---|---:|---:|
| 72 steps, original gate | +36.020441 ppm | +35.917497 ppm |
| 144 steps, 5e-15 local gate | +36.020607 ppm | +35.899087 ppm |

The warm difference is stable across these two numerical profiles.
The absolute warm gain rises by about 7.7 ppm in both arms, while the
absolute cold gain rises by about 31 ppm. Thus matching the expanding
and fixed arms removes most of the numerical-profile shift; do not
treat the absolute gains as equally stable. The local ledger receipts,
the raw strict arms and the comparison script are saved separately.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
p=early_pole_live_scalar_exterior_frw_reference.py
python -W error "$p" "$a" --duration 2.4e-5 --panels 144 --elastic-energy-budget 5e-15
python -W error "$p" "$a" --duration 2.4e-5 --panels 144 --elastic-energy-budget 5e-15 --expands
```

The benchmark uses a prescribed radiation clock and frozen polarization.
It does not solve the cosmological expansion or estimate the late relic
abundance. The small redshift increment is a paired finite-run diagnostic.
