# Cross the hard-aperture and cold-mode refinements in the actual held evolution

Use the same two local interventions on uniform hard depth
1 or depth 4 and scalar support of 8 or 32 cold modes.
Unlike the [initial cross](INDEPENDENT-EARLY-POLE-DAUGHTER-ORIGINAL-INTERVAL-HARD-REFINEMENT-2026-11-27.md),
all four corners are genuine unsplit trajectories. The daughter
interval, exposure, charged-state preparation, both elastic
channels and tolerances are matched as applicable. The refined
32-cold-mode fixture uses the exact nested initial seed and
the published [cold rebin method](early_pole_rebin_held_initial_cold_support.py).

| Hard depth | 8 cold modes: final paired mean | 32 cold modes: final paired mean |
|---:|---:|---:|
| 1 | +5.17245309e-4 | +5.08194351e-4 |
| 4 | +3.42243608e-4 | +3.34406346e-4 |

The mixed finite correction, (cold refinement at depth 4)
minus (cold refinement at depth 1), is +1.21370e-6.
It is small compared with the change from hard refinement
but is resolved in these four numerical runs; the effects
are not identical or assumed separable. At 32 cold modes
the deeper hard aperture reduces every terminal paired
mode contribution, but all 32 remain positive. At each
hard depth the cold-grid change preserves the saved
cumulative sign-crossing bracket. The peak energy defect
across the eight arms is 2.89e-15. This is not a continuum
extrapolation or a cosmological estimate.

Use `early_pole_joined_simultaneous_chart_reference.py` with
`results/early-pole-first-window-hard-depth4-cold32-initializer-2026-11-27.json.gz`
for the new pair at nominal tolerances and
`--record-mode-history`; add `--local-from-paired-archive` for
its local counterpart. Both new raw receipts are retained.
Pass native, then local in each of `--h1q8`, `--h4q8`,
`--h1q32`, and `--h4q32` to
`early_pole_daughter_finite_hard_cold_cross.py`, along with the
earlier initial hard/cold cross JSON as `--static`.
