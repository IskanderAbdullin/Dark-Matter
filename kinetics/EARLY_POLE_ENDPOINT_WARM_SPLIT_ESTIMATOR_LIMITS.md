# The warm split discrepancy is a useful ranking, not a local error bound

At each fixed S4 snapshot, compare full-cell G4 and split-half G4 against **both** accepted Volume32 and Volume64 references, in the same original total-source units. The split estimate is absolute change caused by the split. Even on the snapshots where its rankings correlate strongly with both reference errors, it understates both reference residuals in 6–11 of the 22 complete warm parents. An x=.04 example is q=[0.90829148,1.31278482]: the estimate is 6.95824e-5, but the two reference differences are 8.94897e-5 and 8.86650e-5.

Moreover, using the split result itself makes the source farther from both references in 4–5 parents per snapshot. In the recurrent worst q=[8.28009,11.96750], the residual increases by at least 1.54e-5 of total net source across all four tested snapshots, even taking the smaller of the two reference increases. At x=.05 it rises from [2.3613e-6,1.0831e-6] to [1.9693e-5,2.0971e-5].

This explains why a split ranking should select budgeted candidates, not certify a stopping tolerance, and why retaining an independent integration control matters. Both comparison references are finite; agreement between them does not make these observations rigorous continuum bounds.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-endpoint-warm-source-cells; d=2026-11-28
python -W error early_pole_endpoint_warm_split_estimator_limits.py \
 --states "$p"-x{01,02,04,05}-state-$d.json.gz \
 --references "$w"-x{01,02,04,05}-$d.json
```
