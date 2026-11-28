# Retrospective fixed-budget minimax charts

Use the accepted S4 frozen sources at x=.01, .02, .04 and .05 to minimize
the worst per-epoch sum of absolute net-pressure cell residuals to Volume64.
The integer witness picks exactly one positive order (G4 to G8) in each warm
cell and leaves every cold and charge node unchanged. A separate weighted
Bellman calculation on convex mixtures of the four epochs supplies a lower
bound on the integer optimum. The witness gives the upper bound; the gap is
retained rather than asserting exact integer optimality.

| Scalar node budget | Worst previous .01-trained residual | Worst new feasible residual | Certified lower bound |
|---:|---:|---:|---:|
| 219 | 0.0316515% | 0.0271377% | 0.0270343% |
| 255 | 0.0259624% | 0.0249828% | 0.0249019% |

These are approximately 14.3% and 3.8% reductions in the worst observed
residual, not a uniform improvement at each epoch. Both new charts pay a
higher loss at .01 than the original .01-trained selection. Fresh paired
radial integrations of all eight chart/epoch combinations independently
match the reassembled per-cell values within 3.7e−16 of the full normalized
reference. The paired radial defect is below 2.4e−13. All four snapshots
are used in training, so this is not a new held check, evolved rule,
or a continuum/physical error bound.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results
python -W error early_pole_endpoint_warm_budget_minimax.py \
 --states "$p/early-pole-fine-charge2-scalar4-x"{01,02,04,05}-state-$d.json.gz \
 --reference "$p/early-pole-endpoint-warm-source-cells-x"{01,02,04,05}-$d.json \
 --comparisons "$p/early-pole-warm-frozen-gauss8-source-x"{01,02,04,05}-$d.json \
 --selection "$p/early-pole-warm-budgeted-gauss-probe-$d.json" --direct
```
