# Retained warm moments preserve the reference pressure sign, not every formal probe

For each actual G8, G9 and mixed300 state at x=.04 and .05, fix each warm cell's nonnegative number and mean excess energy measured at the common native mass. Keep the cold cells in their archived constant convention. The calculation allows all positive in-cell measures with these moments. It checks two interpolated energy meshes and includes an analytic curvature bound on the interpolation residual; the reported inner extrema come from explicit two-node mixtures feasible for the original kernel. The finite extrema therefore witness a sign reversal when they straddle zero, independently of the outer numerical enclosure. The two stress kernels are pressure and pseudo-pressure, not transfer functions.

At both checkpoints and at their respective native masses, G9 and mixed300 have more pressure and pseudo-pressure than G8 throughout the admissible sets, and mixed300 remains below G9 for both quantities. At formal mass 0.1 the pressure excesses of G9 and mixed300 over G8 are still strictly positive in these bounded sets. Their pseudo-pressure signs are not determined: at x=.05 the positive measures reach approximately -0.01486% and +0.07286% for G9/G8 pseudo-pressure. The G9 versus mixed300 pressure sign is also not determined at that formal mass. At formal mass 1, both signs are achievable for every pair under either kernel at both checkpoints; for example x=.05 G9/G8 pressure has witnesses at -0.23364% and +0.30377%.

This strengthens the retained-moment interpretation of the *native-mass* endpoint. It does not fix a heavy-mass observable, change masses in the evolution or establish convergence in time or in the momentum grid. These stress masses are mathematical probe parameters, not physical epochs or inputs interchangeable with a cosmological transfer calculation. Outer intervals are analytic-remainder enclosures evaluated in ordinary floating point, not correctly rounded intervals.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-warm
for x in 04 05; do
 python -W error early_pole_endpoint_warm_three_chart_stress_hull.py \
  --native "$p-x$x-state-$d.json.gz" --control "$w-single-g9-upgrade-$d.json" \
  --accepted "$w-mixed300-actual-prefix-x$x-$d.json" \
  --states "$p"-{gauss8warm,gauss9warm,mixed300warm}-x$x-state-$d.json.gz
done
```
