# Outside compensation nearly accounts for the weakened local rebound

Use the same held fine graph and cold seed in four unsplit
calculations. Apply either the exact two-cell aperture-local
daughter swap, the [bare-moment compensation](INDEPENDENT-EARLY-POLE-LOCAL-DAUGHTER-RETUNED-TIME-2026-11-27.md)
on the complementary 18 daughter cells alone, or both changes.
The compensation alone uses the identical affine FD coefficients
as in the combined run while keeping the two aperture cells native.
Its daughter count and bare-energy changes offset those of the
local operation; the combined preparation matches those two native
moments. The comparison enforces that the same compensating
operation appears in both inputs and all terminal gates pass.

| Prepared input | Mean cold-source shift from native |
|---|---:|
| Local only | +5.17245309084e-4 |
| Complement compensation only | -4.28892889953e-4 |
| Local plus complement compensation | +8.85152055818e-5 |

The combined result is within +1.62786e-7 of the sum of the
separate effects, a remainder about 0.038% as large as the
compensation-only effect. In this controlled comparison the
compensation acts as a largely independent negative contribution,
not chiefly through interaction with the two swapped cells.
That statement is finite and conditional on this particular
compensator. It does not identify a unique effect of number,
bare energy, any single elastic process, or the physical spectrum.
The [initial slope check](INDEPENDENT-EARLY-POLE-DAUGHTER-RETUNE-INITIAL-SLOPE-2026-11-27.md)
changes the native-to-local elastic directional difference only
about 0.25%, despite the large finite compensation.

From `kinetics/`, run the simultaneous chart with the same
paired archive and `--compensation-only`. Compare native, local,
compensation-only and combined raw receipts with
`early_pole_daughter_complement_factorial_compare.py`.
