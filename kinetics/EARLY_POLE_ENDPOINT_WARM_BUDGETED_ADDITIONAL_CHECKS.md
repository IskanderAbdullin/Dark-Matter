# Frozen mixed warm orders on two further held snapshots

The 219- and 255-node order maps selected at x=.01 are carried unchanged and directly re-integrated at x=.02 and .04 with both radial meshes. Counts, cold chart and optical fixtures are held fixed at each snapshot. Direct evaluations agree with the independent local assembly to 3.61e−16 of reference net source. The G6/G7/G8 comparators were evaluated separately on the same two states; the old G2–G6 held results reproduce exactly.

| Held x | G6 local L1 | G7 local L1 | G8 local L1 | Mixed 219 | Mixed 255 |
| --- | ---: | ---: | ---: | ---: | ---: |
| .02 | 0.0417017% | 0.0262756% | 0.0243195% | 0.0226477% | 0.0216519% |
| .04 | 0.0480097% | 0.0288385% | 0.0261513% | 0.0316077% | 0.0259624% |

The 255-node mixed rule retains a smaller L1 than uniform G8 in this four-snapshot check, but its margin at .04 is only 0.000189 percentage points. The 219-node rule loses to both G7 and G8 at .04. This is a numerical check across held states on one baseline trajectory, with no re-optimization at .02 or .04 and no guarantee for an independently evolved mixed chart.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results; d=2026-11-28
python -W error early_pole_endpoint_warm_budgeted_additional_checks.py \
 --selection "$p/early-pole-warm-budgeted-gauss-probe-$d.json" \
 --states "$p/early-pole-fine-charge2-scalar4-x"{02,04}-state-$d.json.gz \
 --reference "$p/early-pole-endpoint-warm-source-cells-x"{02,04}-$d.json \
 --comparisons "$p/early-pole-warm-frozen-gauss8-source-x"{02,04}-$d.json
```
