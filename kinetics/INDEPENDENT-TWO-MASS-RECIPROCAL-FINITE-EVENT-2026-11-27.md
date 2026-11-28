# Finite event balance in a reciprocal quadratic medium proxy

The [earlier frozen flavor-mass study](INDEPENDENT-EARLY-POLE-FROZEN-MASS-GATE-2026-11-26.md)
provides a reciprocal, illustrative quadratic gauge functional
for a massless-loop proxy with finite bare-energy denominators.
Use that specified functional for a manufactured finite event.
Let the two hard moments be J=J₀+B n, linear in the active
parent and daughter shell counts, with a held photon moment J_γ.
The active energy proxy is

    E(n)=E_bare·n+A[J_γ Σ_i J_i + ½ Σ_i J_i²],  A=g²/π⁴.

Start with an independently allowed bare two-body event. Represent
the fermion daughter by two positive bracketing shell representatives
and the scalar by its fixed energy. The event direction ν(w) has
daughters summing to one count each and one removed parent.
To conserve the full proxy through a finite signed extent ξ,
choose the weight in the positive interval by

    ν(w)·∇E[n+ξ ν(w)/2]=0.

Because the interaction term is quadratic, the midpoint gradient
is an exact discrete derivative. Keeping the initial frozen
energy-matching weight instead gives positive drift
½ ξ² νᵀHν for the same fixed direction, independent of whether
the update is forward or inverse.

The ordinary and soft-scalar fixtures test both signs and two
extent sizes. With frozen matching, the largest manufactured
drift is 1.224e−8 in fixture energy units. Re-solving the positive
weight closes the maximum residual below 3.4e−19; reversing
from the new state recovers the weight and counts. A one-sided
fixture verifies the condition for rejecting an unavailable
bracket instead of extrapolating. This is energy and number
bookkeeping for isotropic shell counts, not an event-level
vector-momentum or entropy-production check.

The calculation tests one specified proxy. It does not match
finite-mass thermal self-energies, interaction widths, a physical
medium's occupation feedback, or cosmological evolution.
Changing the actual dispersion requires a consistent interaction
functional or another energy ledger derived from the physical model.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_reciprocal_finite_event_reference.py
```
[Full receipt](results/two-mass-reciprocal-finite-event-2026-11-27.json).
