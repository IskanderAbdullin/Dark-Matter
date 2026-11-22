# Independent forward-peaked finite-vector angular audit (2026-11-22)

## Question, scope, and independent method

The fixed-angle Gauss integration of finite-mediator χe/χν full-Pauli
thermal heat was previously found defective (missing a factor of 2 in
the transfer-coordinate Jacobian; failed real CI #44986169496). The
repair is present and HEAD 25de5e55 has passing Scientific CI
[#44997020168](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44997020168).
But a correct coordinate change does **not** prove quadrature convergence.

At one fully on-shell incoming χ+massless-lepton event, compare:
(a) independent adaptive SciPy QUAD over polar CM cos(theta*);
(b) Gauss cosine; and (c) Gauss transfer u. At every polar/32-azimuth
node, reconstruct boosted outgoing χ and lepton energies, evaluate
all four physical Pauli factors and the FULL nonconstant Dirac vector
spin trace with finite spacelike propagator. No vacuum angular
factorization or posthoc symmetrization. This is a **single-event
angular weak heat**, not full 5D integrated rates or Γ_el.

Benchmark: Mχ=0.1, m=0, p=0.9, k=1.2, incoming relative cosine −0.35,
Tχ=0.9, Tlepton=1, G=0.001, mediator mV=0.03,0.1,0.3 in
consistent GeV units. Adaptive polar absolute tolerance 1e−24,
relative 2e−12, maximum 400 subintervals. Its internal error estimates
at the three mV are 2.48e−27, 7.97e−26, and 2.86e−24 respectively;
these are **integration error estimates**, not full physical error bars.

## Exact independent Jacobian control and documentation correction

x=|t|=2p*²(1−c), u=x/(mV²+x), umax=4p*²/(mV²+4p*²).
Therefore

    c(u)=1−mV²u/[2p*²(1−u)]
    |dc/du|=mV²/[2p*²(1−u)²].

The squared propagator becomes (1−u)², canceling the denominator
of the Jacobian. Thus exactly:

    ∫_{-1}^1 dc (mV²/(mV²+2p*²(1−c)))²
      = 2mV²/(mV²+4p*²).

Independent regression verifies this identity for p*²=.001,.1,2 and
mV=.003,.03,.3,3,100 to ~floating precision. The previous prose in
FINITE-MEDIATOR-TRANSFER-ANGLE-AUDIT.md accidentally omitted
the factor of **2** in the denominator; this checkpoint corrects
the prose to match the actually repaired source.

## Measured angular convergence (signed fractional error vs adaptive QUAD)

| mV | Adaptive polar heat | Cosine 24/48/160 | Transfer 24/48/160 |
|---:|---:|---:|---:|
| 0.03 | 1.00334556722e−13 | −2.8746%, +1.9318%, +0.00619% | −0.61784%, +1.05073%, +0.00294% |
| 0.10 | 7.28578668010e−12 | +1.0022%, +0.009872%, ~0 | +0.22371%, +0.004783%, ~0 |
| 0.30 | 2.56142533315e−10 | +0.0001213%, <5e−11%, ~0 | +0.0002168%, <5e−11%, ~0 |

**Negative result:** when mV=0.03, a 48-node transformed Gauss integral
is still **1.05% too high**, and its direct-cosine counterpart is
**1.93% too high**. Even 96 transformed nodes have +0.05096% error.
The error as a function of order is nonmonotonic. The change of
variables helps but is *not* an adaptive convergence certificate.
For mV=.3, 48-point both agree with the independent adaptive reference
at the ~1e−13 level for this specific incoming state.

Reproduction:
- Source: finite_mediator_forward_event_reference.py
- Regression: tests/test_finite_mediator_forward_event_reference.py
- Compact JSON: results/finite-mediator-forward-event-2026-11-22.json
- Command: python -m pytest -q tests/test_finite_mediator_forward_event_reference.py
- Command: python -m kinetics.finite_mediator_forward_event_reference

Next P0: integrate this third adaptive reference over the whole
incoming (p,k,mu) phase space for massive e and massless chiral ν;
independently refine incoming radial, polar, azimuthal quadrature and
thermal tails, vary mV/T and split ν/e bath temperatures; preserve
absolute as well as relative errors when Q nearly cancels. Then
actual mediator/χ1χ2φ gauge content and constraints,
e±γν physical kinetics, nonthermal χ relaxation and fully coupled
Bose production/FRW. No 86-eV ΩDM, BEC, frozen f(q) or CLASS claim.
