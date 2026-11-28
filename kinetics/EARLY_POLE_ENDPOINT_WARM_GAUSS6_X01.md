# First restart from the checked Gauss6 startup

The first G6 continuation from x=.0075 to .01 passed 512/1024 paired time steps on its own unchanged 255/448/488 chart. The saved endpoint binds the exact startup bytes. Its scalar yield is 3.474976525e−7, paired yield drift is 1.25e−7, and absolute continuity residual is at most 1.28e−16.

At the same native S4 mass and the same x=.01 clock, G6 pressure is 0.0801523% lower than G5 and its total number is lower by 1.21e−8 in relative fraction. The independently computed interval of relative pressure differences conditioned on each chart's own stored cell count and mean energy is [−0.000801523682, −0.000801522089]. The similar startup and x=.01 values do not determine subsequent stages or the pressure sensitivity to replaying the first restart on the coarser clock.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_gauss_warm_chain_receipt.py --prefix --order 6 \
 --reports "$p-gauss6warm-start-$d.json.gz" "$p-gauss6warm-x01-$d.json.gz" \
 --states "$p-gauss6warm-start-state-$d.json.gz" "$p-gauss6warm-x01-state-$d.json.gz"
python -W error early_pole_endpoint_gauss_warm_comparison.py "$p-x01-state-$d.json.gz" \
 "$p-gauss6warm-x01-state-$d.json.gz" "$p-warm2-x01-state-$d.json.gz" \
 "$p-warm4-x01-state-$d.json.gz" --x .01
```
