# Evidence for a conservative implicit step in the screened elastic model

The [screened-live explicit control](INDEPENDENT-EARLY-POLE-LIVE-SCREENED-ELASTIC-2026-11-26.md)
needs a small first step and many evaluations on a very short
reference interval. The following test is for that same finite
two-flavor packet kernel with fixed polarization and finite
support. It changes the numerical integration, not the scattering
physics or the archived occupations.

## Implicit elastic substep

The unknowns are admissible Fermi logits. The new moments are
formed by positive Pauli quadrature, not by accepting an
arbitrary unconstrained moment update. Form the variational
gradient at the midpoint of the initial and trial moments, then
correct it along their difference so its dot product with
that difference equals the finite change in total energy. This
is a [discrete-gradient construction](https://web.ma.utexas.edu/users/og/PUBLICATIONS/paper_DisHamSys.pdf)
applied here to the reciprocal field energy. The small interaction
difference is computed directly from the quadratic gauge moments
rather than by subtracting two large total energies.

Use the corrected energy for the reciprocal packet remap and the
final-state entropy variables in a backward step. The individual
fermion counts and this discrete energy close, while the concavity
of the Fermi entropy gives a lower bound on its actual increase:
the step length times the nonnegative packet entropy-production
rate at its final affinities. The analytic Pauli Jacobian used
for damped Newton iterations holds the remap fixed in the
Jacobian only; the actual nonlinear residual recalculates it.
The line search, roundoff energy ledger, and entropy checks are
separate from the approximate Jacobian.

Pure elastic tests on both archived states accept single steps
`8e-8` and `4e-7`, with energy/count closure of order `4e-14`
or smaller and the expected entropy inequality. The step-halving
moment differences are still about `1.6e-5`–`4.1e-5` in absolute
units. A solvable, conservative step need not be an accurate
one for the narrow cold source.

A separate symmetric option also corrects the midpoint gradient
of negative entropy along the finite moment difference. At
convergence this yields the exact finite entropy increase, up
to the nonlinear residual, rather than the backward step's
inequality. Its damped Newton solve uses the converged backward
step as a predictor. This removes observed nonlinear line-search
failures from starting at the archived state but is not evidence
for unconditional solvability or L-stability.

## Split coupled comparison

An explicit cold (and optionally photon-pair) half step on either
side of the implicit elastic step uses the same variational
energy and scalar Hartree terms. The cumulative entropy ledger
records the actual entropy increase of the backward elastic
substep, not just its lower bound. A split with the backward elastic substep is first order despite
the symmetric placement of the outer half steps. The symmetric
elastic option is tested independently rather than assuming an
order of accuracy at a nonsmooth spectral aperture.

At exposure `8e-8` with photon pairs off, the explicit reference
gains are `5.25257e-9` (fine) and `6.55038e-9` (third).
Relative errors on the identical short problem are:

| Elastic step | Split steps | Fine gain error | Third gain error |
| --- | ---: | ---: | ---: |
| Backward | 1 | +10.23% | −5.96% |
| Backward | 4 | +2.83% | −1.45% |
| Backward | 8 | +1.43% | −.69% |
| Symmetric | 4 | −1.59% | −.03% |
| Symmetric | 8 | −.67% | +.11% |

The symmetric rule improves the eight-step comparison, though the
third-archive error is not monotone. Paired-photon tests at four
backward split steps also close the ledgers and keep the gains
positive. At the larger duration `4e-7`, sixteen backward split
steps give `2.67737e-8` and `3.18730e-8`, with combined energy
drifts below `3e-12` and positive sampled increments. This is
a longer numerical demonstration. On the fine archive, the
independent explicit run at the same duration yields
`2.66011e-8` (916 RHS calls); the backward split differs by +.65%.
No comparable explicit long reference has been used on third.

This reference has no phase-space, frozen-field, or FRW
convergence claim. It does not resolve the representation
uncertainty of the initial occupations, or the final relic.
Neither unconditional positive occupations for an arbitrary
step nor a speedup over the explicit solver is asserted.

Reproduce from `kinetics/` (one BLAS thread):

```sh
python early_pole_screened_elastic_discrete_gradient_reference.py fine --dt 4e-7
python early_pole_cold_screened_split_reference.py fine --count 8
python early_pole_cold_screened_split_reference.py fine --count 8 --method midpoint
python early_pole_cold_screened_split_reference.py third --count 4 --pair
python early_pole_cold_screened_split_reference.py fine --duration 4e-7 --count 16
```

Raw: `results/early-pole-screened-discrete-gradient-2026-11-26.json.gz`.
