# Actual G7 continuation through x=.02

The exact G7 x=.01 state is restarted at 1024/2048 steps and independently
verified against the stored raw state and both time traces. At x=.02 it
gives scalar yield 1.02181406994e−5, paired relative yield drift 2.185e−6,
maximum continuity defect 1.680e−13 and maximum number-pair defect
1.025e−16.

At the common native S4 measurement mass G7/G6 count is higher by only
2.67e−9 as a relative fraction, while pressure is higher by 0.0591566%.
The symmetric split between changed counts per cell and changed pressure
per particle is +0.0568698 and +0.00228681 percentage points, respectively.
This keeps the pressure sign and the main contributing warm cell seen
at the first continuation. The conditional pressure interval at fixed
stored count and mean energy does not bound time integration or physics.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_g7_prefix_check.py \
 --native "$p"-{start,x01,x02}-state-$d.json.gz --native-chain "$p-chain-$d.json" \
 --reports6 "$p"-gauss6warm-{start,x01,x02}-$d.json.gz \
 --states6 "$p"-gauss6warm-{start,x01,x02}-state-$d.json.gz \
 --reports7 "$p"-gauss7warm-{start,x01,x02}-$d.json.gz \
 --states7 "$p"-gauss7warm-{start,x01,x02}-state-$d.json.gz
```
