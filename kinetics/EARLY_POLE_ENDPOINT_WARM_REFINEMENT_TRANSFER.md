# An early local error test selects two useful warm cells on later snapshots

At x=.01, rank complete warm cells by the absolute change in their net pressure source when replacing one positive four-node q²-volume rule with the same rule on both halves of that cell. Rank only on that discrepancy, without using Volume64 labels or any later snapshot in the selection. The top two parents are [3.96368,5.72884] and [5.72884,8.28009], identical to the previously identified 219-node chart. Compare against the same two-cell budget assigned by largest G4 net source or largest current population.

| x | Uniform G4 | Uniform G5 | Two cells by split discrepancy | Two largest source cells | Two most populated cells |
| ---: | ---: | ---: | ---: | ---: | ---: |
| .01, selection | 0.06797% | 0.04767% | 0.01612% | 0.06811% | 0.06717% |
| .02, newly evaluated | 0.07190% | 0.04855% | 0.01706% | 0.07203% | 0.07108% |
| .04, newly evaluated | 0.07416% | 0.05593% | 0.02491% | 0.07424% | 0.06616% |
| .05, prior checkpoint | 0.09510% | 0.05502% | 0.03708% | 0.07809% | 0.07744% |

Values are sums of absolute **signed net** cell errors divided by the matching Volume64 total. The reference is paired against Volume32, not assumed exact. Scoring against Volume32 instead of Volume64 also gives the split-discrepancy selection a lower L1 than uniform G5 at every snapshot. In the two new references the whole-band Volume32/64 L1 drift is 1.22e-5 and 1.34e-5; their worst tiny high-q cell fails the separate per-cell gate, so this is a total-normalized error comparison rather than a uniform local accuracy claim. The two newly evaluated snapshots are from the same S4 evolution, not statistically independent simulations. A direct whole-chart evaluation with only the two selected changes exactly reproduces the independent cell replacement predictions on all four snapshots. No guarantee for other grids, states or moments follows. A targeted evolved chart requires its own trajectory and fixed-clock controls.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-endpoint-warm-source-cells; d=2026-11-28
# For a new fixed snapshot, generate its paired reference with
# early_pole_endpoint_warm_source_cells.py on that snapshot.
python -W error -m unittest test_early_pole_endpoint_warm_refinement_transfer
python -W error early_pole_endpoint_warm_refinement_transfer_test.py \
 --states "$p"-x{01,02,04,05}-state-$d.json.gz \
 --references "$w"-x{01,02,04,05}-$d.json
```
