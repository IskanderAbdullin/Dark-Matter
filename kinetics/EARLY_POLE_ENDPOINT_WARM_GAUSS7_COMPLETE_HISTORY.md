# Complete paired G7 warm history under the fixed native cold/charge convention

The 277/448/488-node G7 chart completes the original five paired epochs
through x=.05. Each restart uses the exact preceding refined state.
At the last epoch the scalar yield is 3.71627290723e−4, the paired
yield drift is 3.18010e−5, and the maximum continuity and pair-number
defects are 2.775e−11 and 1.12e−16, respectively.

At the native S4 terminal observation mass, count changes by +8.72e−10
relative to G6. Pressure rises by +0.0615328%; the changed between-cell
counts contribute +0.0591895 percentage points and pressure per particle
+0.00234334. The sign is opposite to the earlier G5→G6 change. The nominal
pressure and the conditional range at fixed cell count and mean energy
both remain positive. An independent 24/48 quadrature of the exponential
reconstructions confirms the terminal stress comparison.

At sampled formal observation masses 0, 1e−6, 1e−5, the native mass,
.01 and .1, the conditional pressure and pseudo-pressure signs are
positive. At masses 1, 3, 10 and 100, finite-atom witnesses within
the conditional count-and-mean sets allow both signs, even though the
exponential reconstruction stays positive. The formal masses are not
physical times. This is neither a convergence claim nor an estimate
of solver error across the full history.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_g7_prefix_check.py \
 --native "$p"-{start,x01,x02,x04,x05}-state-$d.json.gz --native-chain "$p-chain-$d.json" \
 --reports6 "$p"-gauss6warm-{start,x01,x02,x04,x05}-$d.json.gz \
 --states6 "$p"-gauss6warm-{start,x01,x02,x04,x05}-state-$d.json.gz \
 --reports7 "$p"-gauss7warm-{start,x01,x02,x04,x05}-$d.json.gz \
 --states7 "$p"-gauss7warm-{start,x01,x02,x04,x05}-state-$d.json.gz
python -W error early_pole_endpoint_frozen_stress_sweep.py "$p-x05-state-$d.json.gz" \
 "$p"-warm{2,4}-x05-state-$d.json.gz --gauss "$p-gauss7warm-x05-state-$d.json.gz"
python -W error early_pole_endpoint_warm_stress_conditional_hull.py \
 --native "$p-x05-state-$d.json.gz" --gauss "$p"-gauss{6,7}warm-x05-state-$d.json.gz \
 --receipts results/early-pole-frozen-stress-sweep-w124-gauss{6,7}-$d.json --levels 6 7
```
