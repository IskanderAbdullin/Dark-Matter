# Matched mixed255/G7 continuation through x=.02

The fixed mixed255 rule completes the actual 1024/2048 pair from its
verified x=.01 state. Its terminal scalar yield is 1.02181406976e−5 and
its paired yield drift is 2.185e−6. The x=.02 continuity residual stays
below 1.69e−13. Both the original frozen selection hash and the exact
restart state remain bound.

At the common native S4 mass, mixed255 pressure differs from G7 by
−0.0143283%, compared with −0.0146303% at x=.01. The relative total
count difference is only −1.83e−10. The symmetric cell split attributes
−0.0111681 percentage points to the changed counts per cell and
−0.00316021 to changed pressure per particle. This is an evolved numeric
prefix only. It does not show which chart is physically closer to a
continuum solution.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_fixed_budget_contrast.py \
 --control results/early-pole-warm-budgeted-gauss-probe-$d.json \
 --native-states "$p"-{start,x01,x02}-state-$d.json.gz --native-chain "$p-chain-$d.json" \
 --mixed-reports "$p"-mixed255warm-{start,x01,x02}-$d.json.gz \
 --mixed-states "$p"-mixed255warm-{start,x01,x02}-state-$d.json.gz \
 --g7-reports "$p"-gauss7warm-{start,x01,x02}-$d.json.gz \
 --g7-states "$p"-gauss7warm-{start,x01,x02}-state-$d.json.gz --attribution
```
