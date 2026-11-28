# Actual G7 continuation through x=.04

The exact G7 x=.02 state is restarted at 1024/2048 steps and independently
verified against both stored traces and its terminal state. At x=.04,
G7 gives yield 1.73209229325e−4, paired relative yield drift 1.27671e−4,
maximum continuity defect 7.21e−11 and pair-number defect below 9.37e−17.

At the common native measurement mass, G7 pressure is +0.0596997%
relative to G6, while the relative total count contrast is +1.20e−9.
The symmetric terms from cell count and pressure per particle are
+0.0574709 and +0.00222882 percentage points, respectively. The
earlier positive numerical pressure contrast survives this restart;
no time-error bound for pressure follows from the count trace and
no physical extrapolation is implied.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_g7_prefix_check.py \
 --native "$p"-{start,x01,x02,x04}-state-$d.json.gz --native-chain "$p-chain-$d.json" \
 --reports6 "$p"-gauss6warm-{start,x01,x02,x04}-$d.json.gz \
 --states6 "$p"-gauss6warm-{start,x01,x02,x04}-state-$d.json.gz \
 --reports7 "$p"-gauss7warm-{start,x01,x02,x04}-$d.json.gz \
 --states7 "$p"-gauss7warm-{start,x01,x02,x04}-state-$d.json.gz
```
