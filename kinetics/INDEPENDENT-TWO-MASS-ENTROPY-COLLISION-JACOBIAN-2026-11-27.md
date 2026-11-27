# Collision Jacobian and implicit event-space solve on the offset chart

Differentiate the manufactured
[entropy-coordinate collision chart](INDEPENDENT-TWO-MASS-ENTROPY-EXTENT-2026-11-27.md).
For a signed event flux S=F−R, with A=log(F/R), use

    dS = S dlog(F) + R dA.

This retains the signed small difference rather than differentiating
two nearly cancelling large numbers. The event geometry and
stoichiometric weights stay fixed during differentiation. Check the
result at thermal equilibrium and at a displaced parent occupation
by independent fourth-order directional differences.

At detailed balance F=R=M on each event and S=0. Let H be the positive
diagonal negative Hessian of the entropy with respect to the bin
counts. Then the full sum of event Jacobians is

    J = −Σ_j M_j ν_j ν_jᵀ H.

After conjugating by H^(1/2), it is symmetric negative semidefinite.
The deliberately three-event chart has exactly three dissipative
directions. The two individual daughter-plus-parent number
combinations and the split collision energy are left nulls both at
and away from equilibrium. Away from equilibrium the extra
S dlog(F) term generally prevents the same symmetry or sign argument.

Using G_j=√M_j H^(1/2)ν_j, the detailed-balance implicit resolvent
reduces to an inverse on the three event coordinates:
( (I+hGG^mathsf{T})^{-1}
=I-hG(I+hG^mathsf{T}G)^{-1}G^mathsf{T}).
The event-space and dense direct solves agree within 2.23e−14 relative
in the ordinary and near-rest controls, including h=1000 in the
fixture time units. Independent Jacobian differences are below
9.6e−16 relative at equilibrium and 1.49e−11 away. The small
near-rest energy uses the explicit rest-excess term in the
invariant, not a difference of absolute energies.

This is a fixed manufactured chart. No solver performance
claim, nonlinear stability away from equilibrium, physical
source, or cosmological abundance is inferred.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_entropy_collision_jacobian_reference.py
```
[Full receipt](results/two-mass-entropy-collision-jacobian-2026-11-27.json).
