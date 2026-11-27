# Preserve the original local momentum interval under hard refinement

The original local daughter substitution uses two cells covering
the same physical momentum interval that has eight cells on
the fourfold hard mesh. Use the exact recorded continuous
projection in those cells and the original native projection
everywhere else. Match the cold seed and full hard edges within
each paired archive; verify the incoming raw native source.
This evaluates one shared local physical interval, not one
fixed number of cells that changes size with resolution.

The original initial local-minus-native cold rate is -0.01551056,
compared with -0.01167084 after hard refinement. All eight
initial cold-bin shifts remain negative; the total suppression
weakens by about 24.75%. This is a static discretization
sensitivity. It does not by itself test the finite local rebound on
the more resolved hard mesh or imply another choice of physical
momentum interval behaves the same.

Use `early_pole_daughter_original_interval_archive_control.py`
with both compatible fine archive bundles.

## Finite unsplit control at fixed cold resolution

Keep eight cold modes, both elastic channels, the held response, the
same physical daughter interval, and identical solve tolerances and
exposure samples. The only changed discretization is the hard aperture.
The transplant receipt matches the independently checked static
indices, endpoints, and initialized cold rate exactly.

| Compared with native on the same hard grid | Original hard aperture | Fourfold hard aperture |
|---|---:|---:|
| Local daughter: final mean difference | +5.17245309e-4 | +3.42243608e-4 |
| Full continuous archive: final mean difference | -9.75456362e-4 | -1.16341964e-3 |
| Local daughter: cumulative zero bracket | [1.250, 1.458]e-7 | [1.875, 2.083]e-7 |

The positive local rebound survives the harder support test but
weakens by 33.83% and crosses later in the recorded samples. The
whole-archive difference becomes more negative instead. The local
swap is not an additive decomposition of the whole-archive result.
Peak absolute energy closure error among the two refined local
control arms is 1.56e-15. This is a two-grid result, not an extrapolated
hard-mesh limit or a physical relic prediction.

From `kinetics/`, the exact compact refined initializer is
`results/early-pole-first-window-exact-hard-depth4-initializer-2026-11-27.json.gz`.
Use it in `early_pole_joined_simultaneous_chart_reference.py` with
`--local-from-paired-archive` for the refined local arm and without
that switch for its native comparison. Use
`results/early-pole-joined-fine-216-memory-2026-11-27.json.gz` for the
matching original-depth pair. The comparison accepts both full
continuous controls optionally via `--global-continuous` in
original, then refined order.

## Static cross with the cold-output mesh

Separately evaluate the initialized cold RHS on the 8- and
32-mode cold meshes for both hard resolutions. The new cold
cells nest in the original edges, and their archived seed counts
aggregate back to the eight old bins within 8.67e-19. The two
estimated initial scalar poles agree within 2.1e-18.

| Local-minus-native initialized RHS | Eight cold modes | 32 cold modes |
|---|---:|---:|
| Original hard aperture | -0.01551056 | -0.01625060 |
| Fourfold hard aperture | -0.01167084 | -0.01223583 |

All retained modes give negative initialized shifts in the
four fixtures. The coarse-to-fine hard change in the initial
shift is +0.00383972 on the old cold mesh and +0.00401477 on
the finer cold mesh. Their difference is +0.000175049, or
4.56% of the old-cold hard correction. This quantifies a
nonzero initial grid interaction and does not determine the
mixed correction to the finite trajectory. Reproduce with
`early_pole_daughter_initial_hard_cold_cross.py --original ORIGINAL --refined REFINED`
using the two archives named above.
