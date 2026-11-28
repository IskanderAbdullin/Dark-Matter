# Positive warm volume nodes reduce the frozen collision-source bias at equal node count

The [frozen W2 control](INDEPENDENT-EARLY-POLE-ENDPOINT-WARM-PRESSURE-SOURCE-CONTROL-2026-11-28.md) showed a substantial warm-grid effect at the same scalar and charged occupations. Compare warm scalar quadrature rules directly on the **unchanged S4 piecewise-constant input** at x=.01 and x=.05. The charge chart, scalar mass, and clock normalization are held fixed for each epoch. Only complete warm cells are changed; the first has lower edge q=.0109296321, so this comparison excludes the partly straddling q=.01 cell.

The reference repeatedly samples the same parent-cell scalar occupation with 8, 16, 32, and 64 Legendre nodes weighted by q²dq per warm cell. It evaluates the actual collision source and q²/E together at every scalar node. The original one-node and refined geometric rules represent the source as constant within each subcell and weight it by the exact volume-averaged q²/E there. The [positive Gauss2 volume rule](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-POSITIVE-CELL-QUADRATURE-2026-11-28.md) instead samples the source-weight product at its two moment-matched nodes.

| Warm rule | Scalar nodes | Warm pressure source relative to the 64-node reference, x=.01 | x=.05 |
| --- | ---: | ---: | ---: |
| Geometric W1 | 145 | +9.2572% | +9.2213% |
| Geometric W2 | 167 | +2.3975% | +2.4846% |
| Geometric W4 | 211 | +.5524% | +.5410% |
| Positive Gauss2 | 167 | −.03489% | +.08227% |

For the same 167 scalar nodes as W2, the positive rule is much closer to the high-order frozen source at these two epochs. The reference's 32→64 comparison is within 2.63e−6 relative over the recorded gross and signed coefficients, and 4→8 charge-radial tests and event energy projections pass. This does not prove a spatial error bound for arbitrary states or a full evolving pressure calculation. The input shapes were deliberately held fixed. The actual warm tail and spectrum may change under further evolution.

From `kinetics/`, run `early_pole_endpoint_warm_collision_volume_reference.py` on the archived S4 `x01` and `x05` states. The [early](results/early-pole-endpoint-warm-volume-reference-x01-2026-11-28.json) and [late](results/early-pole-endpoint-warm-volume-reference-x05-2026-11-28.json) receipts include the full coefficient and node ladders.
