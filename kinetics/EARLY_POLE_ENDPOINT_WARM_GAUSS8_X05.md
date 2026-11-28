# Complete actual G8 continuation to the planned x=.05 clock

The G8 continuation joins the startup and x=.01, .02, .04 and .05 by exact restart hashes and rechecks the matched G7 history against the same native scalar-mass measurement at each clock. At x=.05 the scalar yield is 3.71627290406e-4; the 1024/2048 time-step drift is 3.18010e-5, the maximum absolute continuity defect 2.7741e-11, and the maximum paired number defect 2.1468e-17.

At x=.05 G8 pressure is 0.0499069% below G7 and total count differs by -8.519e-10 in relative units. The paired-step pressure envelope is [-0.0499149%, -0.0498990%]. Approximately 95.2% of the pressure change is between-cell number redistribution (-0.0475181 percentage points), versus -0.00238888 points within the cells. The original native cell q=[3.96367894,5.72884111] contributes 71.1% of the net pressure difference. This checks the actual whole G8 trajectory to the planned endpoint on one retained spatial chart, not convergence in warm order or a physical matter-power/abundance prediction.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_g7_prefix_check.py --levels 7 8 \
 --native "$p"-{start,x01,x02,x04,x05}-state-$d.json.gz --native-chain "$p-chain-$d.json" \
 --reports6 "$p"-gauss7warm-{start,x01,x02,x04,x05}-$d.json.gz --states6 "$p"-gauss7warm-{start,x01,x02,x04,x05}-state-$d.json.gz \
 --reports7 "$p"-gauss8warm-{start,x01,x02,x04,x05}-$d.json.gz --states7 "$p"-gauss8warm-{start,x01,x02,x04,x05}-state-$d.json.gz
```
