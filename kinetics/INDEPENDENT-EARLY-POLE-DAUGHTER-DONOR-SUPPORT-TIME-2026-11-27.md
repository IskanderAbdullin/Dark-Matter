# The finite local response depends on the shape of equal-moment compensation

In both retuned cases, keep the same two exact aperture-local
daughter cells and restore full retained daughter number and
bare energy to native. In one case all 18 outside cells may
participate in the affine retune; in the other, forbid the
immediately adjacent cell on each side. Verify the
[static receipts](INDEPENDENT-EARLY-POLE-DAUGHTER-EQUAL-MOMENTS-DONORS-2026-11-27.md)
match the evolved initial preparations exactly. Hold hard and
cold grids, numerical tolerances, both elastic channels, model
response and all saved times fixed. All four trajectories pass
the original global checks.

| Preparation | Final mean local-minus-native cold rate |
|---|---:|
| Local without compensation | +5.17245309e-4 |
| Compensate on all outside cells | +8.85152056e-5 |
| Exclude both immediate neighbors from compensation | +3.00531110e-4 |

The last two cases match the same two bare moments yet differ
by +2.12015905e-4 in the finite contrast (about 3.40×).
The initial cold shifts are nearly identical. This demonstrates
finite dependence on where the compensation is put, not a
unique decomposition into density and energy effects or a
prediction for a physical relic.

Use the held simultaneous chart with the fine raw archive and
`--retune-donor-support no_neighbors`; compare with the accepted
native, local and all-outside retuned outputs plus the static
donor receipt using `early_pole_daughter_donor_support_time_compare.py`.
