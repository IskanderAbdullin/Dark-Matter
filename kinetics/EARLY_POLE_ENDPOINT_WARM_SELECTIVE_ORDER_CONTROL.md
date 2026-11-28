# Splitting the selected warm cells works better than merely increasing their order

Keep each original S4 snapshot's constant scalar and charged occupations, scalar mass and collision clock. Compare four-node pressure-source integration on the unrefined cells with three alternatives on the same two predeclared warm parents: (a) five nodes across each entire selected cell, (b) six nodes across each entire selected cell, and (c) separate sets of four nodes on each half. Build each hybrid chart directly; the measured two order variants reproduce the corresponding per-cell substitutions and leave the other integrals unchanged. Compare against the previously accepted 32 and 64 reference-volume rules.

| x | Selective order 5, 213 scalar nodes | Selective order 6, 215 | Split order 4, 219 |
| ---: | ---: | ---: | ---: |
| .01 | 0.047050% | 0.041067% | 0.016117% |
| .02 | 0.048803% | 0.041668% | 0.017059% |
| .04 | 0.056991% | 0.050637% | 0.024913% |
| .05 | 0.071347% | 0.063944% | 0.037083% |

Values are global L1 errors of the signed net per-cell source, scaled by the corresponding Volume64 total. The split is better than either selective order rule for all four snapshots against both reference orders. At .05 uniform six-node integration (255 scalar nodes) gives 0.052998% by the same metric; the cheaper selective sixth-order rule does not match it at that snapshot. The result distinguishes subdivision from putting more nodes in exactly the same cells; it does not establish a dynamical ordering or a continuum configuration.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-endpoint-warm-source-cells; d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_warm_selective_order_control
for x in 01 02 04 05; do
 python -W error early_pole_endpoint_warm_selective_order_control.py \
  "$p"-x$x-state-$d.json.gz "$w"-x$x-$d.json
done
```
