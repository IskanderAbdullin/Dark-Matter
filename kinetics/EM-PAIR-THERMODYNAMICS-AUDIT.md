# Finite e±–photon equilibrium thermodynamics and charge-neutrality gate (2026-11-22)

## Scientific problem, hypothesis and scope

The independent massive χe 5D validation fixes the *conditional scattering
normalisation* but leaves a major physical failure mode: treating the
electron population as an always-massless, inexhaustible thermal bath.
Test the minimum **massive, finite-heat** electromagnetic (EM)
equilibrium reference before attempting χ+SM time evolution.

**Hypothesis:** the μe−=μe+=0 ideal equilibrium e−+e++γ bath tends
from g*ρ,EM=5.5 at T≫me to 2 at T≪me; its finite positive heat
capacity and entropy agree with independent FD/Bessel expressions.
At sufficiently low T, imposing even a tiny net positive background
charge (illustrative η=6×10⁻¹⁰ per photon) requires μe−=−μe+≠0;
a zero-μ e± equilibrium is no longer charge-neutral.

Predeclared failures: Bessel/continuum ρ or P discrepancy >3×10⁻¹¹,
central finite-difference C=∂ρ/∂T relative error >5×10⁻⁷
at T/me=.1,1,10, nonpositive C, or failure to meet neutral-charge
target to 10⁻⁹ relative. Domain .02≤T/me≤100, p/T Gauss 0..48.
This is **equilibrium thermodynamics**, not dynamic annihilation,
instantaneous neutrino recoupling, kinetic detailed balance of
e+e−↔γγ, physical G(T), or a full finite-bath χ collision solver.

## Explicit equations and independent reference

For massive e± together (g=4), x=me/T, q=p/T, ε=√(q²+x²):

    ρ_pairs/T⁴ = g/(2π²) ∫dq q² ε f(ε)
    P_pairs/T⁴ = g/(6π²) ∫dq q⁴/ε f(ε)
    C_pairs,μ=0/T³ = g/(2π²) ∫dq q² ε² f(ε)[1−f(ε)],
    f(ε)=1/[exp(ε)+1].

The last formula means ∂ρ/∂T **at fixed particle mass and μ=0**;
it assumes fast pair chemical equilibrium, not fixed separate electron
and positron numbers. Independent fugacity-Bessel reference:

    P_pairs/T⁴ = gx²/(2π²) Σ_{n≥1} (−1)^(n+1) K₂(nx)/n²,
    ρ_pairs/T⁴ = gx²/(2π²) Σ_{n≥1} (−1)^(n+1)
                 [3K₂(nx)/n²+xK₁(nx)/n].

With photons ργ/T⁴=π²/15 and Pγ=ργ/3, define
g*ρ,EM=30(ρpairs+ργ)/(π²T⁴),
g*s,EM=45(ρpairs+Ppairs+4ργ/3)/(2π²T⁴)
and C_EM/T³=C_pairs/T³+4π²/15.
Neutrinos are **not** put into the electromagnetic heat reservoir;
their temperatures/weak interactions require separate treatment.

For *illustrative* charge neutrality, let η=npositive/nγ=6e−10,
nγ/T³=2ζ(3)/π², and μe−=−μe+=yT (γγ chemical equilibrium).
Solve

    (n_e−−n_e+)/T³ = (2/(2π²))∫dq q²
                       sinh(y)/[cosh(ε)+cosh(y)]
                     = η 2ζ(3)/π².

The stable hyperbolic difference avoids catastrophic subtraction of
two nearly equal FD populations at high T.

## Independently calculated benchmark

Gauss–Legendre p/T quadrature with 256 nodes, cutoff 48 and independent
alternating K₁/K₂ series with 400 terms. Relative Bessel agreement
for ρ and P ~1×10⁻¹⁴ or better in the central cases. With
δT/T=10⁻⁴, independent fixed-m central-difference ∂ρ/∂T
matches integrated C within 1.35×10⁻⁷ at T/me=.1 and much
better near/above T=me.

| T/me | g*ρ,EM | g*s,EM | C_EM/T³ | (n_e−+n_e+)/T³ at μ=0 |
|---|---:|---:|---:|---:|
| 0.05 | 2.00000336683 | 2.00000264208 | 2.63191850594 | 5.13053480e−8 |
| 0.1 | 2.01546141741 | 2.01258974305 | 2.69204866454 | 4.35873611e−4 |
| 0.2 | 2.55570813415 | 2.47803335025 | 3.93063657667 | 2.68411259e−2 |
| 1 | 5.23060273448 | 5.14700573740 | 7.06078240895 | 0.30699373189 |
| 2 | 5.43417936679 | 5.40678447628 | 7.19529259194 | 0.34885213549 |
| 10 | 5.49745687960 | 5.49620334087 | 7.23604196365 | 0.36468185707 |

**Negative massless-bath control:** applying the e± ultrarelativistic
g*ρ,EM=5.5 at T/me=.1 is incorrect: the actual equilibrium
g*ρ,EM≈2.01546. Conversely replacing EM by a photon-only bath
before pair annihilation is invalid near T≈me.

**Negative charge-neutrality control:** if η=6×10⁻¹⁰ is imposed
as an illustrative instantaneous ratio, solving both electron
and positron FD densities gives:

| T/me | μe−/T | n_e+/n_e− |
|---|---:|---:|
| 1 | 5.0819e−10 | 0.99999999905 |
| 0.1 | 3.35315e−7 | 0.99999932938 |
| 0.05 | 0.00284868 | 0.994318844 |
| 0.025 | 13.8376847 | 9.56621e−13 |
| 0.02 | 23.5121362 | 3.78110e−21 |

At T/me=.025, a zero-μ pair bath cannot encode the residual
electron charge; y≈13.84 is required within the model. **This
does not derive a cosmological history of η or e± annihilation**.
Actual η(T), proton content, photon reheating and neutrino decoupling
are separate unknowns.

## Reproduction, limitations and next physical gate

Code: [em_pair_thermodynamics.py](em_pair_thermodynamics.py).
Tests: [test_em_pair_thermodynamics.py](../tests/test_em_pair_thermodynamics.py).
Compact machine-readable [results](results/em-pair-thermodynamics-2026-11-22.json).

    python -m pytest -q tests/test_em_pair_thermodynamics.py
    python -m kinetics.em_pair_thermodynamics --T-over-me 0.1
    python -m kinetics.em_pair_thermodynamics --T-over-me 1.0

The independent closed series is a reference for equilibrium EOS
only. No e+e−↔γγ matrix element/collision source, finite SM neutrino
sector, charge-conserving χ/e reaction dynamics, UV gauge completion,
Friedmann/entropy integration, Bose production or actual freeze-out
has been demonstrated.

Next P0: add *microscopic* pair-creation/annihilation/weak channels
and a charge-neutral finite-temperature evolving e±–γ–ν bath, with
Qdark+QEM+Qν=0, independent number/charge/entropy ledgers,
then couple the full nonlinear χe collision operator (not
phenomenological BGK). Test transfer cutoffs, ν left-doublet
sector and mediator validity above EW restoration.

Background: early-universe quantum collision accounting
(https://arxiv.org/abs/1701.08844); electron-positron annihilation
freeze-out, *Phys. Rev. D* **101**, 063507 (2020),
https://doi.org/10.1103/PhysRevD.101.063507.
No observational 86-eV DM, Ω, BEC or CLASS statement follows.
