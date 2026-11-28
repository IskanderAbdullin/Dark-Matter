# Matched mixed255/G7 continuation through x=.04

The fixed mixed255 rule completes the actual 1024/2048-step pair from its
verified x=.02 state. Its terminal scalar yield is 1.73209229321e−4,
with paired yield drift 1.27671e−4 and maximum continuity residual
7.21e−11. The original frozen rule and all previous checkpoints
remain bound by their original hashes.

At the common native S4 mass, mixed255 pressure is −0.0149062% relative
to G7, with a relative count difference of only −2.40e−11. The
symmetric count and pressure-per-particle contributions are
−0.0117395 and −0.00316670 percentage points, respectively. This is
a numerical comparison of actual prefixes, not a continuum error
estimate or a completed x=.05 history.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_fixed_budget_contrast.py \
 --control results/early-pole-warm-budgeted-gauss-probe-$d.json \
 --native-states "$p"-{start,x01,x02,x04}-state-$d.json.gz --native-chain "$p-chain-$d.json" \
 --mixed-reports "$p"-mixed255warm-{start,x01,x02,x04}-$d.json.gz \
 --mixed-states "$p"-mixed255warm-{start,x01,x02,x04}-state-$d.json.gz \
 --g7-reports "$p"-gauss7warm-{start,x01,x02,x04}-$d.json.gz \
 --g7-states "$p"-gauss7warm-{start,x01,x02,x04}-state-$d.json.gz --attribution
```
