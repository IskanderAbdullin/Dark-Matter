# The completed 219-node targeted history lies between uniform G4 and G5 pressure

Keep the native S4 cold cells and charged D2 geometry. In the 22 complete warm cells, use four positive volume-Gauss nodes except in the two predeclared cells [3.96368,5.72884] and [5.72884,8.28009]. For those, split at the arithmetic midpoint and place four nodes on each half. The genuine independently evolved chart has 219 scalar, 448 daughter and 488 parent nodes. Its metadata records both exact splits and the actual prior uncompressed state hash; a restart with a different split is rejected. The five paired history stages through x=.05 pass the conservation, positivity, mass and timestamp gates.

| x | Target/G4 pressure minus one | Target/G5 pressure minus one |
| ---: | ---: | ---: |
| .01 | +0.078392% | −0.033611% |
| .02 | +0.079455% | −0.033975% |
| .04 | +0.076141% | −0.034340% |
| .05 | +0.078827% | −0.033733% |

All pressure measurements use each epoch's common held W1 mass. At .05 the total yield is 0.000371627290328865 and the target/G4 number difference is +1.72e-9 relative. The two refined cells account for essentially all of the target/G4 pressure increase: their parent-level change is +0.0788243% of G4 pressure, and the signed contribution of the other parents is about +0.00000286%. Under the same count/mean-energy-conditioned pressure interval the target remains strictly between G4 and G5 at all four epochs. This is not a statement that it is more physically accurate than either uniform chart.

G5 has 233 scalar nodes and the target 219. The previous cheaper frozen-source result supports the design choice, but neither lower frozen-source error nor the intermediate evolved value proves convergence. The new independent higher uniform trajectory is a separate test.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
t="$p-targetedwarm"
python -W error -m unittest test_early_pole_endpoint_targeted_warm_history
python -W error early_pole_endpoint_targeted_warm_history.py audit \
 --reports "$t"-{start,x01,x02,x04,x05}-$d.json.gz \
 --states "$t"-{start,x01,x02,x04,x05}-state-$d.json.gz
for x in 01 02 04 05; do
 python -W error early_pole_endpoint_targeted_warm_pressure_compare.py \
  "$p"-x$x-state-$d.json.gz "$t"-x$x-state-$d.json.gz "$p"-gauss{4,5}warm-x$x-state-$d.json.gz
done
```

Use `early_pole_endpoint_targeted_warm_history.py advance --state <output>` for the archived startup convention, then continue through the same recorded stages with `--restart`, `--x-end`, and `--steps`. The chain receipt audits the archived calculations; it does not stand in for a dynamical replay.
