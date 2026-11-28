# The terminal mixed255/G7 instantaneous source ranking survives fixed-node shapes

On the verified terminal states, allow arbitrary nonnegative changes to
the weights of the original warm nodes, while preserving the measured
count and mean energy separately in every native warm cell. Keep each
chart's actual charged state and the chosen high-order reference
unchanged. The affine net response in a cell then spans a sharp interval.
The minimum absolute residual is zero if the interval crosses the
reference, otherwise the nearer endpoint; the maximum is the farther
endpoint. The sum is sharp because cell constraints are independent.

| Reference | Allowed L1 residual, G7 | Allowed L1 residual, mixed255 |
|---|---:|---:|
| Volume96 | [0.0338417%, 0.0338570%] | [0.0250492%, 0.0250637%] |
| Volume64 | [0.0342157%, 0.0342311%] | [0.0258334%, 0.0258479%] |

The ranking holds even for the most favorable admissible G7 weights
and the least favorable mixed255 weights on these fixed node sets.
For Volume96, their separation is at least 8.7779e−5 of the
normalized reference. If a different reference altered each chart's
*normalized residual vector* in L1 by less than 4.38896e−5, the
separation would remain. This condition is not a guarantee that
arbitrary further quadratures or continuum shape reconstructions
have such small effects.

```sh
cd kinetics
d=2026-11-28
python -W error early_pole_endpoint_warm_actual_source_rank_budget.py \
 "results/early-pole-warm-actual-source-reconstruction-$d.json" \
 "results/early-pole-warm-actual-source-moment-hull-$d.json"
```
