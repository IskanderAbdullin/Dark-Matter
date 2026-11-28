# The Gauss6 prefix reaches x=.02 on the checked clock

The next Gauss6 continuation (.01 to .02) passed the 1024/2048 paired time receipt, with the exact previously archived fine endpoint as its restart. At .02 its scalar yield is 1.021814067e−5, relative paired yield drift 2.185e−6, and maximum sampled absolute continuity residual 1.69e−13.

With both G5 and G6 measured at the common native S4 .02 mass, G6 pressure is 0.0804901% lower than G5, and its total scalar count is lower by 3.81e−9 relative. This pressure sign remains in the interval [−0.000804908212, −0.000804893082] allowed at the respective archived per-cell count and mean energy. The actual pressure trend up to .02 is close to the .0075 and .01 values. Later stages remain separate, and the conditioned interval is not a pressure time-convergence test.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_gauss_warm_chain_receipt.py --prefix --order 6 \
  --reports "$p"-gauss6warm-{start,x01,x02}-$d.json.gz \
  --states "$p"-gauss6warm-{start,x01,x02}-state-$d.json.gz
python -W error early_pole_endpoint_gauss_warm_comparison.py "$p-x02-state-$d.json.gz" \
  "$p-gauss6warm-x02-state-$d.json.gz" "$p-warm2-x02-state-$d.json.gz" \
  "$p-warm4-x02-state-$d.json.gz" --x .02
```
