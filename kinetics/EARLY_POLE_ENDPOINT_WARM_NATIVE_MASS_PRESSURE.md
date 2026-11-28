# Native pole-mass remeasurement leaves the warm pressure ordering unchanged

At each of the five accepted clocks, read the hashed G7, G8 and mixed255 node occupations. Reproduce the common-S4-mass measurements exactly, then substitute each chart's own numerical pole mass with no re-evolution. The same count-energy reconstruction is used in both evaluations.

Across 15 measurements, the largest absolute fractional pressure correction is 1.575e-11, while the largest fractional mass shift relative to native S4 is 1.018e-8. Counts are unchanged. Both cross-chart pressure orderings survive at all five measured clocks. Thus choosing the native instead of common mass does not explain the finite differences between these warm layouts. This differs from the cold D2 speed comparison, in which each held chart's substantially different mass affects a conditional speed ordering. Here we do not alter the masses in the collision equation or infer a physical pressure result.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28
python -W error early_pole_endpoint_warm_native_mass_pressure.py --root results \
 --closure "results/early-pole-warm-g8-mixed255-three-chart-$d.json" \
 --mixed-lock "results/early-pole-warm-mixed255-full-prefix-$d.json" \
 --control "results/early-pole-warm-budgeted-gauss-probe-$d.json"
```
