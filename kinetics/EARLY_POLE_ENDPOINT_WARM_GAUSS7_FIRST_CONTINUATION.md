# First matched continuation: G7/G6 pressure contrast stays positive

Both warm rules complete an independently audited 512/1024 continuation from
their own accepted fine x=.0075 states. G7 carries forward the exact restart
bytes and gives scalar yield 3.4749765561e−7, paired yield drift 1.25231e−7,
and maximum absolute continuity defect 9.39e−17 at x=.01.

Evaluated at the same native S4 mass, the G7/G6 pressure change moves from
+0.05908166% at startup to +0.05908024% at x=.01. G7 count at x=.01 changes
by only +8.85e−9 relative to G6. A symmetric cell decomposition at that
epoch gives +0.05678588 percentage points from changed population across
cells and +0.00229436 percentage points from changed pressure per particle.
The largest positive contribution is the native cell q=[3.96368,5.72884].
The pressure sign survives the conditional envelope at fixed count and mean
energy in each cell. The envelope does not bound time error, unobserved
history or physical pressure.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_g7_prefix_check.py \
 --native "$p"-{start,x01}-state-$d.json.gz --native-chain "$p-chain-$d.json" \
 --reports6 "$p"-gauss6warm-{start,x01}-$d.json.gz \
 --states6 "$p"-gauss6warm-{start,x01}-state-$d.json.gz \
 --reports7 "$p"-gauss7warm-{start,x01}-$d.json.gz \
 --states7 "$p"-gauss7warm-{start,x01}-state-$d.json.gz
```
