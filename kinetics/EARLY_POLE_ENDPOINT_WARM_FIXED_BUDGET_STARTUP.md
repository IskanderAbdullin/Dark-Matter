# Actual fixed-budget startup: mixed255 versus G7

The mixed chart keeps the exact order map trained at x=.01 in the earlier
frozen-source test, with 255 scalar nodes and the original 448/488 charge
nodes. This is not the later four-snapshot minimax chart. Its actual
x=.0075 startup finishes the paired 2048/4096 clock with yield
5.7823259832e−8, yield drift 5.11556e−5 and maximum continuity residual
2.92e−16. Both raw artifacts bind the unchanged selection by its hash.

At the common native S4 measurement mass, relative to actual G7 startup,
the mixed chart changes total count by −2.40e−9 and pressure by −0.0146104%.
The latter remains negative throughout the conditional pressure envelope
[−0.000146104346,−0.000146103865] as a relative fraction. This is an
observable difference between numerical trajectories at one endpoint. It
does not rank their physical accuracy or carry the held-source result
over to later evolution.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_fixed_budget_contrast.py \
 --control results/early-pole-warm-budgeted-gauss-probe-$d.json \
 --native-states "$p-start-state-$d.json.gz" --native-chain "$p-chain-$d.json" \
 --mixed-reports "$p-mixed255warm-start-$d.json.gz" \
 --mixed-states "$p-mixed255warm-start-state-$d.json.gz" \
 --g7-reports "$p-gauss7warm-start-$d.json.gz" \
 --g7-states "$p-gauss7warm-start-state-$d.json.gz"
```
