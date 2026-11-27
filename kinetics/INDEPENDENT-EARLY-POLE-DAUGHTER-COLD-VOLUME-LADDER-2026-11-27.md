# Scalar-volume quadrature identifies remaining frozen cold-grid error

The [pointwise hard ladder](INDEPENDENT-EARLY-POLE-DAUGHTER-POINTWISE-HARD-LADDER-2026-11-27.md)
separates the reconstructed daughter occupation error from the direct
comparison of two archived shape prescriptions. At fixed hard depth
1 or 4, now refine the retained cold modes to 8, 16, 32, 64, 128 and
256. Keep each hard moment array exact, and verify that regrouped
cold counts reproduce the old eight-shell seed. On each new mesh
measure the paired initial cold rate on the native frozen packets.

For a separate reference, integrate the signed source contrast in
scalar phase volume (q^3), using the archive's piecewise-constant
cold number density in that volume. Cut at that density's original
edges. Additional cuts occur when a parent or daughter cell's
momentum edge meets the angular boundary (p_d=p_{m parent}pm q)
at the exact fixed scalar energy. The two hard grids give 7 and 17
such extra cuts within the retained cold range. Doubling the scalar
Gauss order and increasing the root scan from 129 to 257 agree
within 8.93e-14 in total rate with these geometric cuts; cuts only at
scalar density boundaries converge substantially more slowly.

| Volume-centroid cold modes | Relative error at hard depth 1 | Relative error at hard depth 4 |
|---:|---:|---:|
| 8 | 4.8983% | 4.8985% |
| 16 | 0.06234% | 0.06251% |
| 32 | 0.30775% | 0.30773% |
| 64 | 0.08182% | 0.08180% |
| 128 | 0.002255% | 0.002232% |
| 256 | 0.002912% | 0.002913% |

The volume-integrated reference is -0.0122530013225 and
-0.0122528987791 on the two hard grids. The nonmonotone centroid
errors are distinct from solver time-tolerance errors and from
the occupation reconstruction error. In particular, 32 modes
have a stable time-tolerance result without necessarily being
more accurate than 16 modes for this separate *initial* pointwise
integral. This does not extrapolate any evolved local contrast,
replace a time-dependent cold closure, or infer a physical abundance.

From `kinetics/`, use `early_pole_daughter_pointwise_cold_ladder.py`
on the original paired 8-mode archive. Use
`early_pole_daughter_direct_cold_volume_quadrature.py` on the original
and the exact depth-4 initializer. Select `--split-angular-edges`
for 8/16 Gauss nodes and 1/2 subpanels, and vary `--root-scan`
from 129 to 257. The plain controls omit the angular switch.
Combine them with `early_pole_daughter_cold_volume_ladder_compare.py`.

## Use identical event-aligned scalar evaluation points on both hard grids

Keep the union of the 24 distinct angular cuts, the retained
eight-shell boundaries, and the native scalar density cuts.
The result has 39 volume intervals and preserves the exact
grouped scalar seed. Evaluate each interval with the same
Gauss rule on both hard meshes; the mesh selection did not
use the measured cold-rate errors.

| Scalar evaluations | Relative error at hard depth 1 | Relative error at hard depth 4 |
|---:|---:|---:|
| 39 (one per interval) | 0.18114% | 0.18114% |
| 78 (two per interval) | 10.752 ppm | 10.752 ppm |
| 156 (four per interval) | 0.570 ppb | 0.565 ppb |

The two-point panel rule is more accurate here than 128 and
256 logarithmic-volume-centroid samples. The count is scalar
integrand evaluations, not edge-search work, wall time or
the number of modes retained in a dynamical solver. The cuts
are located in the frozen starting hard state; moving supports
would have to be handled and checked separately. Run
`early_pole_daughter_shared_cold_angular_sampling.py` with the
original and depth-4 fixtures and the combined volume ladder JSON.
