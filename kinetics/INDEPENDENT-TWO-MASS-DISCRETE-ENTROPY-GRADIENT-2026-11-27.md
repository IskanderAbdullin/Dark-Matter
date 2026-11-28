# An unsplit discrete entropy-gradient step

Use the same manufactured positive chart as the
[coupled backward solve](INDEPENDENT-TWO-MASS-COUPLED-EXTENT-NEWTON-2026-11-27.md).
Instead of evaluating the quantum flux at the endpoint, evaluate
the effective forward and inverse occupations using the entropy
gradient averaged along the count segment from n₀ to n₁.
Its exact thermodynamic property is

    Δ𝒮 = ᾱ · (n₁−n₀) = h Σ_j Ā_j S_j(ᾱ) ≥ 0,

where Ā_j=ν_j·ᾱ. The line-averaged logarithms use their exact
positive-argument integrals and a series for small changes.
A feasible Newton solve enforces all positivity and Pauli bounds
without clipping. The count increments are still in the same
event span, so both daughter populations and split rest-plus-
kinetic energy close by construction.

The averaged-gradient step is symmetric on its smooth solution
branch. On four mesh levels it gives observed final orders
2.0002 (ordinary) and 1.9995 (near-rest) against an independent
unsplit tolerance-controlled reference. At 64 steps the
maximum count errors are 1.33e−5 and 1.02e−5, respectively.
The exact entropy identity is also checked at each accepted step.

A separate severe test starts near the occupation walls and
asks for h=500 in fixture units. The averaged-gradient solve
refuses it in both cases. The ordinary case can use a single
verified backward step, while the near-rest case requires
two backward half-steps. This is controlled fallback, not an
accuracy or existence guarantee for an arbitrarily stiff step.
Second-order behavior is claimed only where the symmetric
solver was used without fallback.

The charts have manufactured capacities and rates. This is
a numerical structure and convergence check, not a physical
production or cosmological transfer calculation.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_discrete_entropy_gradient_reference.py
```
[Full receipt](results/two-mass-discrete-entropy-gradient-2026-11-27.json).
