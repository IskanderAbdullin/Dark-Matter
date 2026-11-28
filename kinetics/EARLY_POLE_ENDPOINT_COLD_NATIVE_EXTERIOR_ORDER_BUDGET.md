# How little unknown exterior response could erase the retained momentum ordering?

Take the verified final sign-targeted fine/third windows and leave
all measured panel counts fixed. Admit an independent unknown exterior
fraction up to the same cap in each chart, normalized to that chart's
total response. Suppose the exterior has q between zero and a formal
upper cap, and is outside the retained [1e−6,1.5e−5] interval.

The adverse comparison puts the fine exterior at zero and the third
exterior at the upper cap. If the retained moment gap at the
independently compatible limits is g>0 and the corresponding upper
probe value is H, the critical fraction is g/(g+H). For RMS use the
second moment before taking its square root. This threshold is sharp
as a limiting witness under the specified support assumptions.

For forward number and the illustrative cap q=.0015, the maximum
fraction per chart is 9.3873e−6 for mean q, 8.4968e−8 for RMS q and
1.8662e−5 for mean speed at the common mass. If q is unbounded, no
positive fraction is guaranteed for the first two; the speed threshold
is still 1.61295e−5 because the limiting speed is one. Both native
meshes reproduce the thresholds. These values are not measurements
of missing response and cannot extend the ordering without independent
bounds on that response.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_endpoint_cold_native_exterior_order_budget.py --root results \
 --accepted results/early-pole-cold-native-moment-sign-refinement-2026-11-28.json
```
