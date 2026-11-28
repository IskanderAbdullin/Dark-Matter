# Hard-aperture refinement moves the native/continuous contrast in the held unsplit chart

Start from accepted native and continuous first-window split
snapshots on the original and fourfold local hard subdivisions.
For each compatible archive, mount the exact initial hard
moments, hard edges and the same eight cold bins. Evolve the
simultaneous charged-logit, held-response chart at the
same tolerances and exposure, not an outer-time split or an
active response run. Cross-check that the new starting cold
source matches the corresponding split archive.

| Hard aperture subdivision | Native mean | Continuous mean | Continuous below native |
|---|---:|---:|---:|
| 1 | 0.06649893898 | 0.06552348262 | 1.466875% |
| 4 | 0.06636534824 | 0.06520192860 | 1.753053% |

The signed native and continuous rate changes are -2009 and
-4907 ppm, respectively. The contrast remains positive but
increases by 0.28618 percentage points with the finer support.
Thus the aperture dependence is not explained entirely by
outer-time interleaving. The largest energy defect is 1.56e-15 and
the original entropy and number gates also pass.

This changes the hard mesh and archive projection, and does not
compare the aperture-local daughter swap on that new hard mesh.
The result is neither a hard-grid convergence proof nor a
physical cold-equivalent production estimate.

Run the simultaneous chart once per projection on a
`branches` bundle containing the accepted native and continuous
depth-four split raw results, and once on the original fine
bundle. Use those two native/continuous pairs and the depth-four
bundle with `early_pole_unsplit_aperture_depth_compare.py`.
