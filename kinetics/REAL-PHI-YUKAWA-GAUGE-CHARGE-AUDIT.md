# Gauge/Yukawa completion obstruction for the original real-boson χ1→χ2+φ model (2026-11-22)

## What is (and is not) being excluded

The project's original [microscopic model](MICROPHYSICS-MODEL.md)
contains **two Dirac fermions χ1,χ2 and ONE REAL scalar φ**,
an exact dark Z2 (χ1 and φ odd, χ2 and SM even), the Yukawa term

    L_Y = − y φ (barχ2 χ1 + barχ1 χ2),

and an **electron-only effective vector coupling to χ1**

    L_contact = G (barχ1 γμ χ1)(bar e γμ e).

Recent χ1–e and χ1–ν scattering pilots replaced this contact
by an illustrative anomaly-free gauged U(1)_(L_e−L_μ) vector Z′,
with χ1 a vectorlike charged Dirac singlet and e_L,e_R,ν_eL,ν_μL
coupled as required by the lepton gauge current. An explicit
2025 example confirms anomaly-free lepton-family differences and a
vectorlike singlet dark fermion
[EPJC 85, 81 (2025)](https://doi.org/10.1140/epjc/s10052-025-13829-3).

**Key question:** Can that candidate preserve the ORIGINAL
renormalizable real-φ Yukawa while giving χ1 a nonzero diagonal
Z′ current but χ2 **no** Z′ current? Under the stated simple
vectorlike charge-eigenstate assumptions, **NO**. This is a
scoped, algebraically exact model-consistency obstruction;
it is NOT a no-go theorem for dark matter, all UV completions,
charged complex scalars, Higgs mixing or the 86-eV benchmark.

## Analytic charge selection proof

Under a continuous U(1)′ gauge transformation,

    χ1 → e^(i q1 α) χ1,
    χ2 → e^(i q2 α) χ2,
    φ  → e^(i qφ α) φ.

The Yukawa monomial φ barχ2 χ1 transforms with phase

    exp[i (qφ − q2 + q1) α].

It is gauge invariant for all α **if and only if**

    qφ − q2 + q1 = 0.

A **single real one-component** scalar cannot carry a
nontrivial continuous U(1) charge: a continuous 1D
real representation maps U(1) into O(1)={±1}; the
connected U(1) maps only to +1. Therefore the original
real φ has qφ=0. The Yukawa selection rule forces

    q1 = q2.

For a **nonzero** χ1 vector coupling q1 gχ ≠0,
χ2 necessarily has the *same* nonzero gauge charge.
Its diagonal vector current

    gχ q2 Z′μ (barχ2 γμ χ2)

is gauge invariant and is ALSO Z2-even
(χ2 is even, so its bilinear remains even). It
cannot be omitted from the minimal gauged model by
the original dark Z2. In particular, selecting q1≠0,
q2=0 and a single real φ gives Yukawa residual q1≠0:
**the displayed scalar decay operator is forbidden**
under the attempted U(1)′ completion.

This is independent of all numerical 5D scattering
quadratures, mediator masses, Tν/Tγ assumptions and
the dark scalar 86-eV mass.

## Explicit finite arithmetic and escape routes

The reproducible checker
[gauge_charge_yukawa_audit.py](gauge_charge_yukawa_audit.py)
uses exact rational charges and the U(1)′ phase
selection rule, not floating tolerance. Four
[regressions](../tests/test_gauge_charge_yukawa_audit.py)
cover the analytic statement and representative
q1,q2∈{−3,…,+3} controls:

| χ1 charge | χ2 charge | scalar | Yukawa allowed? | χ2 vector coupling? |
|---:|---:|---|---|---|
| +1 | +1 | real, qφ=0 | yes | **yes** |
| +1 | 0 | real, qφ=0 | **no** | no |
| +1 | 0 | complex, qφ=−1 | yes | no; but changes field content |
| 0 | 0 | real, qφ=0 | yes | absent for both |

The complex charged escape has *at least two real*
scalar degrees of freedom and requires new gauge-breaking,
charge, mass, scalar and dark-Z2 stability accounting.
A symmetry-breaking/mass-mixed fermion scenario may
also yield nontrivial effective couplings, but is
not the simple diagonal-charge model tested here.
Adding a new singlet, charged Higgs, higher-dimensional
operator or different parity changes the microscopic
model and invalidates blindly reusing its original
partial width, degrees of freedom and closed kinetic
charges.

The exact Z2 is **not** itself inconsistent:
χ1 odd × χ2 even × φ odd gives Z2-even Yukawa.
The obstruction is simultaneous U(1)′ **charge**
selection and an incorrectly truncated χ1-only
gauge-mediator current.

## Consequence for research priority

The direct χ1–SM heat and all fitted Γ_el
proxies are **incomplete as a UV-completed three-species
model** unless χ2–SM scattering is also included
for q1=q2≠0, with its own distribution, Pauli
factors and energy exchange. An abundant daughter
χ2 after χ1→χ2+φ could exchange substantial energy
with e±/ν; its rates cannot be assumed negligible
without a quantitative comparison. The actual
χ1↔χ2+φ Bose production, mother depletion, daughter
replenishment and the finite e±γν heat reservoir
must then evolve **together**, conserving physical
charges and total energy, before computing late fφ(q).

A proper next step is to specify a single complete
gauge+Z2 Lagrangian with consistent χ1/χ2/φ charges
and mediator mass generation, derive **both**
χ1ℓ→χ1ℓ and χ2ℓ→χ2ℓ (and any allowed inelastic)
collision operators, cross-check their normalization,
and test physical thermal/gauge/astrophysical bounds.
If χ2 cannot scatter, state a different UV model
and rederive the Yukawa/gauge interactions; do not
continue with mutually incompatible assumptions.

This result is **not** a terminal falsification of the
86-eV Bose dark matter hypothesis. No physical
Γ_el, equilibrium, ΩDM, macroscopic BEC, frozen f(q)
or CLASS-ready observable is established.
Reproduction: python -m pytest -q tests/test_gauge_charge_yukawa_audit.py
