# Time width and hard aperture support interact in the finite joined comparison

Cross 16 and 32 macro steps with the original aperture and the
fourfold subdivision of the two original cold-active hard cells
per charged species. Within a fixed aperture and preparation, the
hard moments, hard edges, cold seed and initial cold-mode rates
match exactly across time widths. In all arms the cold mesh,
active response rule, first-window duration (3.3333333333333335e-7),
elastic quadrature and existing local energy gates are unchanged.

| Macro steps × aperture subdivision | Native mean | Continuous mean | Continuous below native |
|---|---:|---:|---:|
| 16 × 1 | 0.06694243427 | 0.06528878623 | 2.47025% |
| 16 × 4 | 0.06680149350 | 0.06507936701 | 2.57798% |
| 32 × 1 | 0.06672255625 | 0.06540677346 | 1.97202% |
| 32 × 4 | 0.06658244639 | 0.06513950571 | 2.16715% |

At fixed time steps, hard refinement lowers both means. Its
effect is almost independent of time width in the native
preparation: the signed interaction is +8.31e-7. In contrast,
for the continuous preparation the finer time width increases the
mean by +1.17987e-4 on the original hard mesh and by only
+6.01387e-5 on the fourfold mesh. Their difference, -5.78485e-5,
is nearly half of the original hard mesh's time-refinement
effect. Consequently a time correction measured on the original
hard support cannot be carried over unchanged to the finer support.
The largest per-step energy defect is below 4.3e-15 in the eight
accepted legs, which does not imply convergence in either variable.

Use `early_pole_cold_screened_split_reference.py fine` with
`--duration 3.3333333333333335e-7 --count 16` or `--count 32`,
`--pair --screening-update full --aperture-grid fine --snapshots`,
`--compton-strength 1 --same-soft-strength 1 --conversion-strength 1`,
`--elastic-energy-budget 5e-15`, aperture subdivision 1 or 4,
and either `--projection steps` or `--projection cell-conservative`.
Compare accepted raw outputs in the row order above using
`early_pole_first_window_time_aperture_cross_compare.py`.
This local numerical test does not establish a nested support limit
or a cosmological rate.
