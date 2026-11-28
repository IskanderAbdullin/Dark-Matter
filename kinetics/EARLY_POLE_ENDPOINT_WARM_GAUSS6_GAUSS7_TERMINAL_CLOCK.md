# Recovered G7 coarse terminal pressure, compared on the two original clocks

Recover the actual G7 1024-step continuation from its verified x=.04
state. Compare the two G7 pressure states with the corresponding archived
G6 pressure states at the single native x=.05 observation mass. Each
order retains its own exact original x=.04 history. The independent
trace gates and all endpoint/restart SHA links reproduce the paired runs.

G7/G6 total pressure differs by +0.000615323821 on 1024 steps and
+0.000615328123 on 2048 steps, as relative fractions. Moving between
these clocks changes the contrast by only +4.30196e−9. Neither the
sign nor the magnitude of the terminal G7/G6 contrast is explained by
this finite comparison of the last clock segment. This is not a bound
on an earlier epoch or on the complete history.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4
python -W error early_pole_endpoint_warm_g7_first_clock_check.py --terminal \
 --native "$p-x05-state-$d.json.gz" --prefix "results/early-pole-warm-gauss7-full-prefix-$d.json" \
 --raw-reports "$p"-gauss{6,7}warm-x05-$d.json.gz --starts "$p"-gauss{6,7}warm-x04-state-$d.json.gz \
 --coarse "$p"-gauss{6,7}warm-coarse-x05-state-$d.json.gz --fine "$p"-gauss{6,7}warm-x05-state-$d.json.gz \
 --coarse-replays "$p"-gauss{6,7}warm-pressure-clock-$d.json
```
