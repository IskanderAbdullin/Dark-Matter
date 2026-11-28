# How much independent cold shape can reverse the G7/G6 stress sign?

Use the exact G6/G7 x=.02 states. Warm cell count and mean reference
energy remain fixed. In each cold cell the two charts may now assign an
equal fraction of their respective counts to independent shapes, with
the remaining population having a common normalized shape. No cold
particle is allowed outside its original cell. No cold energy constraint
is imposed in this released family.

For a monotone kernel with cold cell bracket [A,B], counts n6 and n7,
the difference between the worst shared and worst independent lower
bounds equals min(n6,n7)(B−A). The common shape chooses the adverse
endpoint according to the sign of n7−n6; the independent comparison
puts G7 at the lower and G6 at the upper endpoint. At fixed release
fraction, combine these two contributions. Bounds on the warm terms
from analytic curvature and monotone dependence on the observation
mass control the full adaptive mass intervals.

* At the actual common reference mass 0.0002635394, pressure and
  pseudo-pressure remain positive even with fully independent cold
  shapes under the stated support constraints.
* At all formal masses in [0,.1], both signs survive 3.5% cold release
  in every cell. The 221 certified intervals tile the range.
* At formal mass 0.000002, an explicit compatible 3.8% release reverses
  pseudo-pressure (difference −3.88202e−14) while pressure stays positive.
  The cold witness puts point masses strictly inside cells; its warm
  states use the original positive nodes. Against a shared witness
  that stays positive, this places the uniform release threshold for
  both signs between 3.5% and 3.8%, without claiming the exact threshold.

This is a fixed-state, conditional uncertainty test. No formal mass is
assigned to a time or physical scenario.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_warm_g7_cold_shape_budget.py \
 --native "$p-x02-state-$d.json.gz" \
 --states "$p"-gauss{6,7}warm-x02-state-$d.json.gz \
 --prefix results/early-pole-warm-gauss7-prefix-x02-$d.json \
 --cover results/early-pole-warm-gauss7-continuous-mass-sign-x02-$d.json
```
