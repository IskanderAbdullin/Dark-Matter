# Prepared daughter controls in the unsplit held-response system

Use the same archived frozen elastic graphs, common initial cold
seed and native hard support. Change only the daughter preparation,
either in the two initially cold-active hard cells or in the other
18 retained daughter cells. The unsplit wrapper checks that the split
input snapshots have accepted energy and entropy closure, that the
unsplit model uses their exact cold seed and support, and that their
initial cold mode rates are reproduced before integrating. Use
DOP853 in the admissible charged-logit chart at 3e-9 relative
tolerance over the first 3.3333333333333335e-7 interval.

| Preparation | Mean cold source | Shift from native |
|---|---:|---:|
| Native | 0.06649893898 | 0 |
| Local daughter only | 0.06701618429 | +0.00051724531 |
| Outside daughter only | 0.06572727839 | -0.00077166059 |

All three pass the original conservation gates; the largest absolute
energy closure defect is 8.9e-16. The outside reduction retains its sign
and nearly its size compared with the earlier active-screening 16-step
split test (-0.00074183628). The strictly aperture-local effect does not:
it is positive in this unsplit control rather than the small negative
effect (-0.00001002991) at 16 active split steps. Thus the earlier
claim that the local finite difference was tiny applies to that particular
split experiment, not to the unsplit result. Do not attribute the
disagreement to a single factor from these data; time splitting
and response treatment both differed. Matched held 16- and 64-step
local controls are required for that comparison. The positive local
and negative outside values also cannot be added to infer a whole-species
perturbation without a separate test.

From `kinetics/`, use the accepted raw local and outside split runs as
inputs to the following. They supply the initial snapshot only;
their time histories are not used by the unsplit solver.
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
b=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
p=early_pole_joined_simultaneous_chart_reference.py
python -W error "$p" "$b" > chart-native.json
python -W error "$p" "$b" --prepared-input split-local-daughter.json > chart-local.json
python -W error "$p" "$b" --prepared-input split-outside-daughter.json > chart-outside.json
python -W error early_pole_prepared_unsplit_compare_reference.py chart-native.json chart-local.json chart-outside.json
```

This is a finite, held-response numerical control with manufactured
shape inputs. It does not establish physical production, time convergence
of the split local effect, or robustness as the hard aperture changes.
