# Fourfold elastic radial order reduces the last observed quadrature shift

Hold the original hard aperture, eight cold bins, 16 macro steps,
duration and active screening fixed. Refine the radial order of
both elastic packet graphs together, at held transfer order one
and angle order two. In each preparation the initial hard moments,
cold seed and cold-mode rates remain exact.

| Elastic radial order | Native mean | Continuous mean | Continuous below native |
|---|---:|---:|---:|
| 1 | 0.06694243427 | 0.06528878623 | 2.47025% |
| 2 | 0.06700051953 | 0.06543588591 | 2.33526% |
| 4 | 0.06698923649 | 0.06542742274 | 2.33144% |

From order two to four the native mean moves by -168.4 ppm and the
continuous mean by -129.3 ppm, smaller than the earlier shifts
of +867.7 and +2253.1 ppm. The contrast moves by only -0.00382
percentage points at this last refinement. It is a narrower
radial-only diagnostic, not a universal convergence or error
statement: transfer, angle, hard support and time discretization
have not been refined here, and neither sequence is monotone.

Both new trajectories pass the previous closure gates. The
largest local elastic defect in the two order-four calculations
is 2.23e-15, under the 5e-15 budget. Reproduce with
`early_pole_joined_elastic_packet_refinement_reference.py` at
`--unlike-quad 4 1 2 --same-quad 4 1 2` and the original split
settings. Use the comparison script's optional `--radial4` pair
alongside the default and radial-two raw results.
