# Reduced extents for overlapping, nonparallel energy brackets

Extend the [parallel-bracket solve](INDEPENDENT-TWO-MASS-PARALLEL-MEDIUM-EXTENT-2026-11-27.md)
to two signed events with different affine directions. Their fermion
destinations use the low–middle and middle–high shell pairs, so a
common scalar quadratic energy coordinate is no longer available.
Use the same class of prescribed reciprocal quadratic energies, not
a fitted thermal sector.

For fixed proposed extents ξ, solve the weights w by the
[certified box contraction](INDEPENDENT-TWO-MASS-SHARED-MEDIUM-WEIGHT-CERTIFICATE-2026-11-27.md).
With v_j=a_j+w_j d_j and the common midpoint gradient g_m, the
individual energy residual is F_j=v_j·g_m. Its linearizations are

    G_jk = ∂F_j/∂w_k = δ_jk d_j·g_m + ξ_k v_jᵀH d_k/2
    X_jk = ∂F_j/∂ξ_k = v_jᵀH v_k/2.

Because G is the positive slope matrix times I−T_w at the fixed
point, a strictly contracting root map T also gives invertibility.
The weight tangent is −G⁻¹X. The endpoint count tangent includes
both the direct event directions and the ξ-weighted bracket-direction
derivatives. Apply this derivative to the line-averaged entropy
variables and the affine entropy factors, then perform a two-variable
extent solve. This does not certify admissibility for unspecified
extents; the box is audited for every evaluated trial in the fixture.

Across the three opposing signed steps, the nested solve agrees
with an independent four-unknown solve (explicit forward and reverse
statistical factors) to 3.5e−18 for extents and 3.5e−17 for weights.
The reduced residual Jacobian agrees with fourth-order finite
differences within 4.2e−13 absolute. Endpoint energy tangent errors
are below 4.1e−20 in the fixed-extent checks; both daughter counts
and occupancy bounds are checked. Integrating over 0.5 prescribed
time units gives order 2.0010 on the last two meshes against DOP853
of the independently stepped instantaneous matching equations.
The maximum endpoint count error falls from 1.08e−4 at four steps to
1.66e−6 at 32. Global energy defects are at most 1.53e−20 and
per-step entropy identities close within 2.46e−18.

This verifies nonparallel bracket coupling and differentiability on
a manufactured chart. It does not establish convergence under mesh,
momentum, or microscopic collision-operator refinement.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_nonparallel_medium_extent_reference.py
```
[Raw receipt](results/two-mass-nonparallel-medium-extent-2026-11-27.json).
