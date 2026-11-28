# Two selectively refined warm parents improve the frozen source with fewer nodes

At the same fixed S4 scalar and charge occupations, mass and clock, start from G4 and split each of two previously identified influential complete warm parents at its arithmetic midpoint: [3.96367894345,5.72884111204] and [5.72884111204,8.28009053085]. Use four positive q²-volume points on each half. This uses 219 scalar nodes compared with 211 for uniform G4 and 233 for uniform G5. The untouched-cell collision ledgers are unchanged exactly on both tested endpoints. Pair the same two radial orders and use the saved Volume64 total as reference.

| x | Uniform G4 L1 per-cell discrepancy | Uniform G5 discrepancy | Targeted discrepancy |
| ---: | ---: | ---: | ---: |
| .01 | 0.06797% | 0.04767% | 0.01612% |
| .05 | 0.09510% | 0.05502% | 0.03708% |

The L1 discrepancy sums absolute **signed net** errors across native cells and divides by reference total net pressure source. At .05 the targeted global signed error is −0.007615%, versus +0.05452% for uniform G5 and −0.09352% for uniform G4. L1 is the stronger comparison because opposing residuals cancel in the signed total. The two target cells were identified from prior frozen/evolved diagnostics, so the numbers do not imply a general grid-selection rule or independent predictive validation. Volume64 is a finite reference, not a rigorous continuum bound. No targeted trajectory has been evolved here.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
for i in 01 05; do
 python -W error early_pole_endpoint_warm_targeted_source_refinement.py \
   "$p-x$i-state-$d.json.gz" "results/early-pole-endpoint-warm-source-cells-x$i-$d.json"
done
```
