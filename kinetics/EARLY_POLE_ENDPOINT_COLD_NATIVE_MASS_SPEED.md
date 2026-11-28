# The native mass change weakens the conditional cold speed ordering

At x=.05, reuse the certified exact fine and third D2 panels on q=[1e-6,1.5e-5]. The third chart's numerical scalar mass is 0.48955% smaller. The common-mass comparison isolates the change in momentum placement. The native-mass comparison evaluates each chart with its own `v(q)=q/hypot(q,m)`. Both bounds allow arbitrary placement in each final panel and reproduce on both quadrature meshes.

| Weight profile | Third minus fine, common mass | Third minus fine, native masses |
| --- | ---: | ---: |
| Forward number | [-3.87714e-4, -1.61297e-5] | [-3.67027e-4, +5.30916e-6] |
| Inverse number | [-3.87866e-4, -1.60492e-5] | [-3.67124e-4, +5.44661e-6] |
| Forward centered excess | [-4.14003e-4, -3.61279e-5] | [-3.91924e-4, -1.32690e-5] |

With native masses, the forward and inverse number profiles no longer determine a sign. The negative forward-excess speed difference still holds. This does not change the earlier common-mass momentum ordering. It is a conditional numerical finding about this retained band only.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28
python -W error early_pole_endpoint_cold_native_mass_speed.py --root results \
  --archive "results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json" \
  --lock "results/early-pole-cold-native-moment-sign-refinement-$d.json"
```
