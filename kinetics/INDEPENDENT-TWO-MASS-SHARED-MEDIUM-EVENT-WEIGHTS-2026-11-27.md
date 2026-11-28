# Individually energy-balanced events need not add conservatively in a shared medium

Extend the [quadratic medium single-event control](INDEPENDENT-TWO-MASS-RECIPROCAL-FINITE-EVENT-2026-11-27.md)
to two manufactured competing reactions sharing the parent and
the fermion daughter's bracketing shells. The two scalar
destinations differ, and one event goes forward while the other
goes in reverse. The medium and all occupations here are test inputs.

For a quadratic energy E with Hessian H, two increments v and w
satisfy the cross identity

    E(n+v+w)−E(n)
      = [E(n+v)−E(n)] + [E(n+w)−E(n)] + vᵀHw.

Even when each increment independently balances its own
medium change from the initial state, their sum can miss the
mixed interaction term. A sequentially recomputed update is
a different scheme; this test compares two simultaneous
updates using isolated versus shared weights.

Solve the two signed extents and the two positive energy
weights using the **same** full midpoint state. Each event
direction is orthogonal to its full midpoint energy gradient.
Use the line-averaged entropy variables of the actual common
before/after counts for both quantum affinities and rates.
Then the total entropy change is the sum of the two
nonnegative affinity-times-rate products, integrated by the
same discrete identity used in the one-event chart.

The three solved fixture steps close both number ledgers,
the full energy to 2.72e−20 in absolute units, and the entropy
identity to 5.3e−18. Independently energy-matched steps
from the old state, combined with those same extents, drift
by −2.99e−7 at the largest tested step. The mixed quadratic
term reproduces that drift. This is an algebraic interaction
effect, not evidence of a physical dissipative sign.
No kinetic convergence, matched thermal dispersion,
actual reaction rate or cosmological prediction is claimed.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_shared_medium_event_weights_reference.py
```
[Full receipt](results/two-mass-shared-medium-event-weights-2026-11-27.json).
