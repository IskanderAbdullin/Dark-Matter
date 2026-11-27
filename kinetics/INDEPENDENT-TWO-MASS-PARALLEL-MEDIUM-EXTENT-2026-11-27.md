# Eliminate parallel energy weights before solving signed event extents

The general [weight box and contraction test](INDEPENDENT-TWO-MASS-SHARED-MEDIUM-WEIGHT-CERTIFICATE-2026-11-27.md)
allows different energy-deposition directions. This narrower case has
an exact reduction: two opposing events share an affine direction d
for their daughter energy interpolation. Their energy variation
depends on only one aggregate weight coordinate. The Hessian H here
is the prescribed constant [quadratic medium](INDEPENDENT-TWO-MASS-SHARED-MEDIUM-EVENT-WEIGHTS-2026-11-27.md).

For proposed signed extents ξ, set D=Σ ξ_j a_j and z=Σ ξ_j w_j, so
the count increment is D+zd. Energy conservation requires

    (dᵀHd/2) z² + [dᵀ(g₀+HD)] z + Dᵀg₀ + DᵀHD/2 = 0.

Use the cancellation-free root next to zero on the positive-slope
branch. Recover w_j=−a_j·g_mid/(d·g_mid) from the shared midpoint
gradient and check every weight and occupation separately. The
quadratic would not by itself certify the individual brackets.

Differentiating the quadratic at fixed initial state gives
dz/dξ_j=−a_j·g_end/(d·g_end). This is the endpoint energy tangent.
Differentiate the recovered weights along it and propagate both
the change in the energy weights and the change in the line-averaged
entropy variables into the 2×2 extent Jacobian. In this way the
extent solve does not have to iterate the two weights as extra
unknowns. The step is still coupled through the common state.

The reduced solve reproduces the published extents and weights
to 5.3e−17 absolute with maximum finite-difference Jacobian
defect 3.1e−13. The energy tangent defect is below 1.4e−20.
Across 4, 8, 16, and 32 steps a separate high-order integration of
the instantaneous energy-matched flow gives observed order 2.001
on the last pair, full trajectory energy defects under 1.0e−20,
and step entropy identity defects under 1.7e−18. This is an
actual convergence result for the prescribed five-state test,
not evidence for kinetic convergence or a physical thermal sector.

Replay from kinetics with the previously published two-event
JSON in results:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_parallel_medium_extent_reference.py
```
[Raw receipt](results/two-mass-parallel-medium-extent-2026-11-27.json).
