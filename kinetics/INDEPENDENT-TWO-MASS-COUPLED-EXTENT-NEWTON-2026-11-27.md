# A coupled implicit step solved without splitting event directions

The [analytic collision Jacobian](INDEPENDENT-TWO-MASS-ENTROPY-COLLISION-JACOBIAN-2026-11-27.md)
allows a small nonlinear solve on the same manufactured chart.
Let V have the fixed energy-preserving stoichiometric vectors as
columns. Solve backward Euler for all event extents together:

    ξ = h S(n₀+Vξ),    n₁=n₀+Vξ.

Newton's matrix is I−h(∂S/∂n)V. Use a feasibility-preserving line
search and stop only on a residual check. Continuation increases h
from a smaller numerical trial to the requested value with the
**same initial state**; these are not time-integration substeps.
The solver refuses a step if it cannot make a residual-reducing
feasible update.

Number and split rest-plus-kinetic energy are linear invariants
of V. By entropy concavity, a feasible solution also obeys

    Δ𝒮 ≥ ∇𝒮(n₁)·(n₁−n₀) = h Σ_j A_j(n₁) S_j(n₁) ≥ 0.

Here the bound is checked independently against the entropy
difference. It does not require a theorem that Newton finds a
unique solution of arbitrary large steps.

Test an ordinary and a near-rest chart, each with a forward,
inverse or almost-empty daughter start, and h=.08, 1, 50 in
fixture time units. The nonlinear step uses three event coordinates
instead of 21 occupation coordinates. All 18 solves retain strict
positive populations and fermion occupancies below one. Largest
event-coordinate residual is 2.33e−11. An independent bounded
occupation-space solve at h=.08 agrees within 1.14e−12 in
maximum count. The normalized energy defect is below 7.7e−18.
These are numerical fixtures, not a claim about arbitrary initial
states, runtime savings on a physical collision operator, or
actual cosmological dynamics.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_coupled_extent_newton_reference.py
```
[Full receipt](results/two-mass-coupled-extent-newton-2026-11-27.json).
