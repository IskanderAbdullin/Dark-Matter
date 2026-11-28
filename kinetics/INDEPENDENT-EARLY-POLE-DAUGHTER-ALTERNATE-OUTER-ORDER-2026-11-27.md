# Centering the slow operator moves the fixed-inner interleaving bias

Compare two symmetric operator orders at 16 macro intervals and
the same 64 individual elastic solves at fixed inner width. The
original order is slow half, elastic (four substeps), slow half.
The alternative is elastic half (two substeps), slow full, elastic
half (two substeps). Keep the held dielectric, native/local inputs,
cold support, exposure and per-elastic step budget fixed and
include the original 64×1 comparison.

| Profile | Paired final mean | Paired first common block gain | Slow RHS calls per arm |
|---|---:|---:|---:|
| Original slow-elastic-slow | +2.12757801e-4 | -1.77275718e-10 | 532 |
| New elastic-slow-elastic | +4.80464418e-4 | -9.71141188e-11 | 308 |
| Fine 64×1 reference | +3.84570917e-4 | -1.25734360e-10 | 1828 |

The two coarse orders fall on opposite sides of the fine
fixed-inner comparison. The simple average of the two coarse
paired means is about +3.46611110e-4, coinciding within 2.45e-10
with the measured original-order 32×2 mean. This is one
observed finite relationship, not a theorem about order
averaging. The fine reference still includes a finite elastic
step. Lower slow RHS calls are a count, not a wall-time estimate
or universal scheme preference. The original gates and each
individual elastic gate pass for the compared paths.

The option only supports held-response experiments. Run
`early_pole_daughter_alternate_outer_order.py` for the native/local
pair with 16 macro intervals, then validate opposite orders
with `early_pole_daughter_alternate_outer_order_compare.py`.
