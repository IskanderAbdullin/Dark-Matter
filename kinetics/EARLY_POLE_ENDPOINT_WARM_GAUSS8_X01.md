# Actual paired G8 first continuation

The exact G8 x=.0075 state advances through paired 512/1024 orders to
x=.01. Its scalar yield is 3.47497653189e−7, the relative paired yield
drift is 1.25231e−7 and the maximum absolute continuity defect stays
below 7.95e−17. All saved clock and restart bytes remain bound.

With the same native S4 mass, G8 count is lower than G7 by 6.98e−9
as a relative fraction, and total pressure is lower by 0.0477456%.
The conditional pressure interval remains negative. The signed
symmetric per-cell pressure terms are −0.0454175 and −0.00232810
percentage points from the changed counts and changed pressure per
particle, respectively. The sign is unchanged from the G8 startup.
This is not yet a later G8 history or a continuum error bound.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_g7_prefix_check.py --levels 7 8 \
 --native "$p"-{start,x01}-state-$d.json.gz --native-chain "$p-chain-$d.json" \
 --reports6 "$p"-gauss7warm-{start,x01}-$d.json.gz --states6 "$p"-gauss7warm-{start,x01}-state-$d.json.gz \
 --reports7 "$p"-gauss8warm-{start,x01}-$d.json.gz --states7 "$p"-gauss8warm-{start,x01}-state-$d.json.gz
```
