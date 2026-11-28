# The fixed mixed255 pressure contrast remains negative after the first restart

Both the fixed 255-node selection and uniform G7 complete the first paired
continuation, at 512/1024 steps from their own exact x=.0075 startup
states. The mixed chart has yield 3.4749765519e−7, paired drift 1.25231e−7,
and continuity defect below 2.39e−16 at x=.01.

At the same native S4 measurement mass, mixed255/G7 pressure changes from
−0.0146104% at startup to −0.0146303% at x=.01; count at x=.01 is
−1.21e−9 as a relative fraction. The x=.01 pressure contrasts for the
changed cell counts and the changed pressure per particle are
−0.0114571 and −0.00317314 percentage points, respectively. The largest
contribution is from q=[5.72884,8.28009]. Each pressure sign remains
negative across the conditional count/mean-energy interval. These are
numeric comparisons, not a ranking of physical accuracy or a full-history
bound, and the trajectory uses the previously frozen .01-trained rule,
not the newer retrospective minimax map.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_fixed_budget_contrast.py \
 --control results/early-pole-warm-budgeted-gauss-probe-$d.json \
 --native-states "$p"-{start,x01}-state-$d.json.gz --native-chain "$p-chain-$d.json" \
 --mixed-reports "$p"-mixed255warm-{start,x01}-$d.json.gz \
 --mixed-states "$p"-mixed255warm-{start,x01}-state-$d.json.gz \
 --g7-reports "$p"-gauss7warm-{start,x01}-$d.json.gz \
 --g7-states "$p"-gauss7warm-{start,x01}-state-$d.json.gz --attribution
```
