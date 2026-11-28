# Warm counts and energy fix the held-mass sign, but not heavier stress probes

Freeze the actual x=.05 G3 and G4 endpoints. In every refined warm cell preserve exactly its archived number and mean excess energy at the common W1 measurement mass. Keep all unrefined cells in their archived constant convention. The previous frozen stress sweep used one positive exponential reconstruction. Here we test all nonnegative within-cell measures with those same two retained moments.

Work in a normalized energy coordinate inside each cell. For either pressure kernel, extrema of its piecewise-linear interpolant use one or two sampled energy points. Each such extremum is a constructive atomic witness for the original (noninterpolated) kernel. To enclose the full continuous class, add the analytic linear-interpolation remainder to each bound. That remainder comes from the second derivative of the kernel using q² as an intermediate variable; use the exact simpler formulas at the reference mass and at zero mass. The code checks independent 64/256 interpolation meshes and both saved nodal and exponential witnesses. Bounds use ordinary floating-point arithmetic, not directed rounding.

At the **reference mass**, G4/G3 minus one for pressure is in [−0.125398%, −0.125382%]; the negative sign is robust to the warm reconstruction. Both pressure and pseudo-pressure remain negative on the tested formal mass grid through 0.1. At **formal M=1**, constructive finite witnesses have both signs: pressure reaches approximately −0.4193% and +0.1171%, while pseudo-pressure reaches −1.0468% and +0.6692%. At **formal M=10**, the pressure witnesses reach roughly −2.8556% and +2.3941%, and pseudo-pressure −12.670% and +13.217%. Each range is for that observable separately, not a single common perturbation optimal for both. An outer interval crossing zero would not establish ambiguity by itself; the actual positive atomic witnesses do.

The heavier test masses are probe parameters applied to a frozen spectrum, not time steps, physical late-universe predictions or a CMB/CLASS result. A chosen exponential fit can yield a definite heavy-mass sign, but the stored warm count and mean energy alone cannot enforce it.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_warm_stress_conditional_hull
python -W error early_pole_endpoint_warm_stress_conditional_hull.py \
 --native "$p"-x05-state-$d.json.gz --gauss "$p"-gauss{3,4}warm-x05-state-$d.json.gz \
 --receipts results/early-pole-frozen-stress-sweep-w124-gauss{3,4}-$d.json
```
