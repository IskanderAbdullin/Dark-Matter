# Spurious conserved neutrino endpoint populations in energy-matched Z′↔ν+antiν event lattices (2026-11-22)

## Scientific question

The first positivity/entropy-preserving nonthermal
ν/antiν/Z′ event graph conserves **all three
expected reduced physical invariants** — energy,
C=Nν+Nantiν+2NV and
L=Nν−Nantiν — but that does **not**
guarantee a unique physical long-time
thermal equilibrium. On a finite discretized
reaction network, modes not connected to
any collision preserve their occupations
artificially. Do the original integer daughter
grid and the later half-offset angular-overlap
grid have exactly the same conserved subspace?

This is a **finite-state reaction topology audit**,
not a new physical ν-equilibration measurement,
no collision amplitude fit, no full SM/FRW
source, and no 86-eV DM conclusion.

## Exact microscopic stoichiometric condition

Each positive on-shell quantum reaction
Z′_j ↔ ν_i + antiν_l
has stoichiometric column

    S_r = (−e_V,j, +e_ν,i, +e_antiν,l).

A linear invariant z is conserved by **all**
possible reaction histories iff

    S^T z = 0.

The three expected physical invariant vectors are

    C=(2 for V, 1 for ν, 1 for antiν),
    L=(0 for V, +1 for ν, −1 for antiν),
    E=(E_V,j, kν,i, kantiν,l).

Each exact lattice event conserves
E_V,j=kν,i+kantiν,l, so all three vectors
are in the left nullspace. Any additional
null direction, if present, represents
a **nonphysical extra invariant** of the
selected discrete network (or a real
model-selection restriction that must be
explicitly justified). Positivity and a
local H theorem by themselves cannot
eliminate such a frozen degree of freedom.

For each lattice, construct integer
species×reaction incidence matrix S.
Its rank is evaluated independently by
finite-field Gaussian elimination
modulo two different primes 1000003
and 1000033, and by singular-value
decomposition over the reals. The three
physical null-vector residuals are
explicitly tested. No floating threshold
alone defines the numerical rank.

## Reproducible strict comparison

Same physical mV=.4, Emax=4,
energy steps h=.4,.2,.1, full angular
kinematic support and positive event
weights. Compare:

1. **Original integer neutrino grid**
   k_i=i h, i=1,...,Emax/h.
2. **Half-offset daughter grid** with exact
   angular-cell intersection,
   k_i=(i+1/2)h, i=0,...,Emax/h−1.

Both have the same number of vector
and neutrino/antineutrino discrete species
for each stated h. They differ only in
the energy collocation and on-shell
reaction connectivity, not in the
physical C/L/E definitions.

| h | States (V+ν+antiν) | Integer events | Integer rank / nullity | Half-offset events | Half-offset rank / nullity |
|---:|---:|---:|---:|---:|---:|
| .4 | 29 | 45 | 24 / **5** | 54 | 26 / **3** |
| .2 | 58 | 189 | 53 / **5** | 207 | 55 / **3** |
| .1 | 116 | 774 | 111 / **5** | 808 | 113 / **3** |

All numbers are independently confirmed
by the two distinct modular primes and
SVD. For h=.4,.2, an additional local
exact-rational SymPy rank calculation
agrees: integer ranks 24,53 and
shifted ranks 26,55, respectively.
SymPy is **not** a required project
runtime dependency; the committed
implementation recomputes modular ranks
without it.

## Explicit origin of two artificial conservation laws

The integer lattice includes the last
neutrino state kν=Emax and last
antineutrino state kantiν=Emax.
For a vector E_V≤Emax, any decay needs
two **positive-energy** daughters with
kν+kantiν=E_V. Hence neither highest
daughter state can appear in any
physical Z′→ν+antiν pair.
The corresponding two incidence
rows are identically zero, and the
two endpoint occupation numbers
are **individually conserved for all
times and every timestep**, regardless
of how accurate or implicit the
numerical solver is.

This explains why the original
integer-grid reaction system has
**five**, not three, conserved linear
combinations. These two additional
null modes are artifacts of the
finite interval/lattice representation,
not new neutrino symmetries.
A long-time solution could maximize
entropy only subject to **five**
constraints, and therefore need not
equal the physically expected
maximum-entropy state determined by
C/L/E alone. Merely reducing dt cannot
remove this structural bias.

On the tested half-offset grid, the
highest daughter energy is Emax−h/2,
which can pair with a lowest positive
daughter energy h/2 at vector
E_V=Emax. The exact angular support
allows this pair at mV=.4 and these
resolutions; no neutrino/antineutrino
row is isolated. The measured
left nullity is **exactly three**,
spanned by C/L/E, on these finite
networks.

**Important limitation:** nullity3
at three synthetic resolutions does
not establish global continuum
ergodicity, correct finite-volume
thermal rates, or absence of
additional approximate slow modes
in a different physical mV/Emax
regime. Quantum Bose/Pauli and
time-dependent mass/SM channels
could alter the process network.
Likewise the split ν and antiν
distributions need not actually
be described by one common
temperature absent weak scattering.

## Decision and reproducibility

The half-offset, angular-overlap
lattice is **topologically preferable**
for this restricted channel: it
removes two rigorously identified
artificial conservation laws of the
original integer grid, independently
of the numerical source error
improvements. This is a scoped
mathematical correction and a useful
negative test to retain in the archive.

- [Executable rank audit](vector_neutrino_event_invariants.py)
- [Focused tests](../tests/test_vector_neutrino_event_invariants.py)
- [Compact integer ranks](results/vector-neutrino-event-invariants-2026-11-22.json)
- [Permanent stage archive](../archive/stages/2026/11/2026-11-22-vector-neutrino-invariants.md)

    python -m pytest -q tests/test_vector_neutrino_event_invariants.py
    python -m kinetics.vector_neutrino_event_invariants

Next verify real CI and perform a
**global constrained equilibrium/long-time
time integration** on the shifted network
against independent Lagrange-multiplier
FD/BE stationary predictions, including
distinct conserved lepton asymmetry.
Then solve remaining physical cutoff
and dt errors, weak ν–e±γ collisions,
gauge-consistent Z′/χ1/χ2 and real
φ Bose FRW. Do not promote any
synthetic Γ0=1 lattice result
to χ Γelastic, 86-eV ΩDM,
BEC, frozen f(q) or CLASS-ready.
