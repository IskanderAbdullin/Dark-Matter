# The actual G8 startup reverses the preceding pressure-order increment

The independently paired 299/448/488-node G8 startup reaches x=.0075
at 2048/4096 steps with yield 5.78232591471e−8, relative paired yield
drift 5.11556e−5 and maximum continuity defect 1.760e−16.
At the same native S4 observation mass, compared with G7, the relative
count change is −1.42e−8, but pressure changes by −0.0477491%.
The population and pressure-per-particle contributions in the
symmetric cell split are −0.0454220 and −0.00232704 percentage points.
The pressure interval at fixed stored cell count and mean energy
retains that negative sign.

G7/G6 had a positive startup increment. The opposite G8/G7
increment rules out using a monotone order trend at this stage.
These are paired startup results only, not a complete G8 trajectory,
monotone convergence statement, or physical pressure.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_g7_prefix_check.py --levels 7 8 \
 --native "$p-start-state-$d.json.gz" --native-chain "$p-chain-$d.json" \
 --reports6 "$p-gauss7warm-start-$d.json.gz" --states6 "$p-gauss7warm-start-state-$d.json.gz" \
 --reports7 "$p-gauss8warm-start-$d.json.gz" --states7 "$p-gauss8warm-start-state-$d.json.gz"
```
