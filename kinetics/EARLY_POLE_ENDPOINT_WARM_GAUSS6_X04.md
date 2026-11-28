# The Gauss6 prefix reaches x=.04

The Gauss6 continuation on x=.02–.04 passes the matched 1024/2048 trace check and links the declared SHA of the earlier fine x=.02 state. The x=.04 scalar yield is 0.0001732092291. The paired relative yield drift is 1.27671e−4 and the maximum sampled absolute continuity residual is 7.21e−11, within the existing 3e−10 gate.

At the native S4 mass for this epoch, G6 pressure is 0.0801953% lower than G5 and G6 total count is lower by 1.43e−9 in relative fraction. Conditioning on the two orders' separately archived per-cell count and mean gives the relative pressure contrast interval [−0.000802001553, −0.000801904026]. This checks the integrated path to x=.04, not the final stage or earlier pressure errors.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_gauss_warm_chain_receipt.py --prefix --order 6 \
  --reports "$p"-gauss6warm-{start,x01,x02,x04}-$d.json.gz \
  --states "$p"-gauss6warm-{start,x01,x02,x04}-state-$d.json.gz
python -W error early_pole_endpoint_gauss_warm_comparison.py "$p-x04-state-$d.json.gz" \
  "$p-gauss6warm-x04-state-$d.json.gz" "$p-warm2-x04-state-$d.json.gz" \
  "$p-warm4-x04-state-$d.json.gz" --x .04
```
