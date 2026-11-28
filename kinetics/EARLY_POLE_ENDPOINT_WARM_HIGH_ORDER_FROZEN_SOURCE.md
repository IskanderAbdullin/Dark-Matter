# Higher positive orders reverse the held G8 warm-source bias

At the two exact, previously certified S4 states x=.01 and .05, hold the charged distributions, scalar number in every original cell, pole mass and clock fixed. Extend the positive q²-volume quadrature to G9, G10 and G12, checking the positive moments and two radial meshes. The G8 aggregate and cellwise absolute net residuals reproduce the prior archive exactly on both held states.

At x=.05 the signed aggregate source residual relative to Volume64 is -0.0275876% (G8), +0.00787851% (G9), +0.0108379% (G10), and +0.0148923% (G12). The same G8-to-G9 sign reversal occurs at x=.01. The leading 8→9 cell is [3.96368,5.72884] on both snapshots; the leading later change moves to [5.72884,8.28009]. The normalized absolute local net difference between G10 and G12 is still 0.00966845% and 0.0122142% at the two clocks. A small aggregate difference alone would conceal this cellwise drift.

Neither Volume64 nor G12 is asserted to be exact. This is an instantaneous source diagnostic on held states, not an evolved G9 history, a causal time attribution, or a convergence proof.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results; d=2026-11-28
python -W error early_pole_endpoint_warm_high_order_frozen_source.py \
 --states "$p/early-pole-fine-charge2-scalar4-x"{01,05}"-state-$d.json.gz" \
 --references "$p/early-pole-endpoint-warm-source-cells-x"{01,05}"-$d.json" \
 --priors "$p/early-pole-warm-frozen-gauss8-source-x"{01,05}"-$d.json"
```
