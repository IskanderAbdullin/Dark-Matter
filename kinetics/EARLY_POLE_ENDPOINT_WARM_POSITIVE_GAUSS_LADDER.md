# Compare positive volume Gaussian orders without evolving the state

At fine S4 x=.01 and .05 hold the occupation and optical/clock parameters fixed. Apply two, three, and four Gaussian nodes for the q²dq phase-space volume only to the 22 complete warm cells; leave the other scalar cells and all charged cells unchanged. On u=(q−a)/(b−a), we solve the symmetric generalized eigenproblem for normalized monomial Gram and multiplication matrices. Both node and weight positivity, exact enclosed support, and polynomial moments through degree 2n−1 are verified. The generalized order-two rule reproduces the prior analytic two-node per-cell rates.

| x | Order | Scalar nodes | Signed total pressure-source error | Sum of absolute cell errors / net reference |
| --- | ---: | ---: | ---: | ---: |
| .01 | 2 | 167 | −0.03489% | 0.10782% |
| .01 | 3 | 189 | +0.05749% | 0.09489% |
| .01 | 4 | 211 | −0.06646% | 0.06797% |
| .05 | 2 | 167 | +0.08227% | 0.20487% |
| .05 | 3 | 189 | +0.04983% | 0.11469% |
| .05 | 4 | 211 | −0.09352% | 0.09510% |

Errors are evaluated against the fixed-occupation 64-node comparison, not a certified bound on the continuum. The absolute cell sums fall at both epochs; the signed total is not monotonic and can change sign. In particular Gauss4 has a *larger* absolute signed error at x=.05 than Gauss3 even while it reduces the absolute distribution of cell errors. Do not select a production order from a single cancellation-prone signed total. The 32/64 reference passed the aggregate check, not a uniform cell-relative gate at high q. These comparisons say nothing by themselves about evolved pressure or a physical spectrum.

From `kinetics/` with the endpoint dependencies:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_warm_positive_gauss_ladder
for i in 01 05; do
 python -W error early_pole_endpoint_warm_positive_gauss_ladder.py \
  results/early-pole-fine-charge2-scalar4-x${i}-state-2026-11-28.json.gz \
  results/early-pole-endpoint-warm-source-cells-x${i}-2026-11-28.json
done
```
