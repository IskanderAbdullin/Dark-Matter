# Complete paired history for the original fixed mixed255 selection

The 255/448/488-node chart completes all five paired epochs from the
original frozen .01-trained order map. The map is **not** the later
retrospective minimax selection. Actual terminal yield is
3.7162729072865e−4; paired yield drift is 3.18010e−5. The continuity
and paired-number defects at the last epoch stay below 2.775e−11 and
5.99e−17, respectively. All restart digests match.

Relative to uniform G7 at the same native terminal observation mass,
count changes by only +1.58e−11 and pressure by −0.0162258%.
Changed counts per cell and pressure per particle contribute
−0.0129338 and −0.00329201 percentage points. The conditional
count-and-energy pressure contrast remains negative across its
[−0.0001623380,−0.0001621784] range in relative units. This verifies an
actual numerical tradeoff with 22 fewer scalar nodes. It neither ranks
physical accuracy nor assumes the held-source ordering is a convergence
statement.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_fixed_budget_contrast.py \
 --control results/early-pole-warm-budgeted-gauss-probe-$d.json \
 --native-states "$p"-{start,x01,x02,x04,x05}-state-$d.json.gz --native-chain "$p-chain-$d.json" \
 --mixed-reports "$p"-mixed255warm-{start,x01,x02,x04,x05}-$d.json.gz \
 --mixed-states "$p"-mixed255warm-{start,x01,x02,x04,x05}-state-$d.json.gz \
 --g7-reports "$p"-gauss7warm-{start,x01,x02,x04,x05}-$d.json.gz \
 --g7-states "$p"-gauss7warm-{start,x01,x02,x04,x05}-state-$d.json.gz --attribution
```
