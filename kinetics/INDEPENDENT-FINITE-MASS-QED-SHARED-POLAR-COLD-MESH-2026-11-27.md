# Separate cold scalar cell weights from bare-kernel node sampling

Use the [shared polar momentum scan](INDEPENDENT-FINITE-MASS-QED-SHARED-POLAR-MOMENTUM-SCAN-2026-11-27.md)
to test a weighted average over q=1e−6 to 1e−5. The test weight
is the frozen archived scalar occupation times q²/E_s(q).
The kernel is the convention-dependent bare finite part at radius
equal to the aperture half-width. Neither weight nor kernel is
being added to a physical rate in this comparison.

The independent fully split integral separates the native scalar
occupation edges and the parent/daughter aperture crossings.
A 17-to-31 node comparison agrees in the weighted average within
1.2e−13 absolute. Its values are −502.09023 (fine) and −417.74262
(third).

Then sample the kernel at the retained q-bin convention, the
center in q³, using 8, 16, 32, 64, and 128 log-spaced cells.
Compare:

1. Point-sample the scalar occupation at the same center and
   multiply by the accurately integrated q²/E_s measure.
2. Keep one kernel evaluation per cell but integrate the true
   archived scalar weight across every native occupation edge inside
   each cell.

On cold32, point sampling overestimates the integrated scalar weight
by 7.45% and 8.80%. Its errors in the normalized kernel average
are −0.274% and +0.221%; the exact per-cell scalar weights reduce
these to −0.0494% and +0.0305%. The exact weights preserve the total
scalar weight independently of the cold output mesh. Nodal sampling
does not show monotonic convergence in the normalization; even
cold128 leaves a third-grid weight error of about −1.78%.
This is a comparison of quadratures for a frozen diagnostic,
not evidence that the time-evolved cold mode operator uses or
inherits these point-sampling errors.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_shared_polar_cold_mesh_reference.py
```
[Raw mesh and normalization results](results/finite-mass-qed-shared-polar-cold-mesh-2026-11-27.json).
