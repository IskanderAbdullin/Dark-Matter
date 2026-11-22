# Finite ν thermodynamic reservoir plus mode-resolved Z′ Pauli/Bose 1↔2: conservative time evolution (2026-11-22)

## Hypothesis, numerical acceptance and scope

A fast microscopic ννbar↔Z′ collision and an analytic closed-system equilibrium do NOT by themselves constitute a time-dependent solution. The next gate is to integrate the actual momentum-dependent quantum source with a **finite** ν energy and particle-number reservoir (not an external fixed-temperature thermostat), and check charge, energy, entropy, time and grid convergence. Use a clearly declared **constrained thermal ν closure**, not an invented fully nonthermal ν Boltzmann operator.

Reduced model: Z′ is an unpolarized massive 3-polarization vector with positive momentum-resolved Bose occupation bV(p,t), two active massless ν flavors plus antineutrinos gν=4, initially FD T0 and zero μ, and bV(p,0)=0. Only Z′↔ννbar 1↔2 operates; all e±γ, χ1/χ2, SM weak number changes, expansion, gauge Higgs and thermal masses are absent. A fast but **UNMODELLED** number- and energy-conserving ν self-redistribution is assumed to keep the ν sector exactly FD with evolving Tν(t), common effective particle-number μν(t) for ν and antineutrinos. This is a specified moment-closure assumption to be tested later, NOT the physical ν thermalization history.

Acceptance: physically nonnegative bV and FD occupations; each substep joint relative Nν+2NV and Eν+EV conserved <2e−12, total entropy nondecreasing; late time at τ=Γν,0 t=3000 reaches independent constrained equilibrium T/T0 and μ/T within 1e−5; independent 96-node CM angular FD blocking at nonzero μ agrees with closed analytic logarithm to 2e−13; explicit failure for an oversized step, and first-order timestep bias recorded (not hidden).

## Exact physical kinetic source and finite ν projection

For each vector momentum p with E=√(p²+mV²), massless ν and antineutrino energies ε±=(E±p cosθ*)/2 and common y=μν/Tν, the ν angle-averaged quantum Pauli coefficient is

    A(p,Tν,y) = 1−<fν(ε+) + fνbar(ε−)>,

    A=1−(2Tν/p) ln{[1+exp(y−(E−p)/(2Tν))]/
                          [1+exp(y−(E+p)/(2Tν))]}.

Now the exact mode source is

    dbV/dτ = (mV/E) A [bBE((E−2μν)/Tν)−bV],
    τ=Gammaν,0 t,   Gammaν,0=g′²mV/(12π).

**This is microphysical Z′↔ννbar** in a prescribed instantaneous FD state. It is NOT a χ elastic collision Γel, arbitrary BGK or a true nonthermal ν source.

For one time substep Δτ at temporarily fixed Tν,y, the mode integrates analytically as

    bV,new = bBE,old + [bV,old−bBE,old]
                         exp[−Δτ(mV/E)A_old].

This is positivity preserving for A≥0 in the tested chemical branch. The updated vector number NV and energy EV are integrated with physical gV p²dp/(2π²). To conserve the genuine collision charge and total energy **without thermal reset**, re-solve Tν,new and y_new from

    Nν(Tnew,ynew)=Ntotal,initial−2NV,new,
    ρν(Tnew,ynew)=Etotal,initial−EV,new.

If either target becomes nonpositive for an overlarge step, abort rather than fabricate a negative ν density or clamp the result. Reprojecting the finite ν bath this way conserves its particle number and energy constraints, but assumes ν kinetic equilibration instantly. Entropy of ν FD and all vector modes is recomputed; if total entropy drops in a substep, abort.

## Numerical results and negative time controls

Use mV/T0=.03, ν g=4 and Z′ g=3, fixed physical vector momentum grid Gauss n96, cut pmax=35T0, neutrino FD integration Gauss n96, cut qν=50, no expansion. Total τ=Γ0t is dimensionless, so **no physical g′/H** is required to compare time convergence.

| τtotal | substeps | Tν/T0 | μν/Tν | NV/T0³ | S_total/S_initial |
|---:|---:|---:|---:|---:|---:|
| 100 | 20 | 1.24573956 | −1.06102130 | 0.05244343 | 1.07001435 |
| 100 | 40 | 1.24469363 | −1.05706073 | 0.05227502 | 1.06998460 |
| 100 | 100 | 1.24405791 | −1.05465133 | 0.05217245 | 1.06996596 |
| 500 | 100 | 1.24892891 | −1.13561898 | 0.06054534 | 1.07284065 |
| 1000 | 100 | 1.24657871 | −1.13627058 | 0.06130994 | 1.07289637 |
| 2000 | 200 | 1.24603652 | −1.13593823 | 0.06142937 | 1.07289861 |
| 3000 | 300 | 1.24600366 | −1.13590414 | 0.06143497 | 1.07289863 |
| 5000 | 500 | 1.24600016 | −1.13589988 | 0.06143550 | 1.07289863 |

The independent **thermodynamic constrained solution**
from the previous separate radial/fugacity audit is
T_final/T0=1.24600008376,
μν/T=−1.13589977975,
S_total/S_initial=1.07289862970.
The dynamic τ=5000 solution approaches it to
|ΔT/T0|≈7.22×10⁻⁸ and |Δ(μ/T)|≈9.88×10⁻⁸.
This is a genuinely independent equilibrium endpoint
control for the time evolution.

Every accepted step conserves the joint N and E
to floating-point root-solver precision; across
τ=3000, n96, steps300 the largest relative
N residual is ≈1.40×10⁻¹⁴ and E ≈1.97×10⁻¹⁴,
with strictly nondecreasing total entropy for
the tested trajectory. There is **no
fixed external ν temperature** or injected energy.

**Negative timestep convergence result:** at τ=100,
the Lie-style frozen-bath exponential update yields
Tν/T0=1.24573957 (20 steps),
1.24469364 (40), 1.24416424 (80),
with a roughly first-order factor 1.98 in
successive differences. Exact conservation and
positive entropy therefore DO NOT imply
time convergence. More decisively, τ=3000 in
only 20 substeps is rejected because a single
stiff step would transfer more particle number
to vector states than the finite neutrino reservoir
has available: the algorithm raises
“step too large or lost bath” rather than
clamping. This is a numerical splitting
failure, **not** a physical exclusion.

**Grid/cutoff controls:** at τ=3000/steps300,
vector n64→96→128 and cut 35T0 give
Tnu/T0≈1.2460036595, 1.2460036591,
1.2460036591; cut25→35 shifts
Tν/T0 by ≈4.32×10⁻⁷ and cut35→45
by <4×10⁻¹⁰. The small cutoff/grid
differences in this smooth thermal ν closure
do not validate arbitrary nonthermal ν
or high-momentum dark production.

At nonzero μν/T=−1.1359, an independent
96-node two-body angular integral of
1−fν(ε1)−fνbar(ε2) matches the analytical
A(p,T,y) to relative 2e−13 at four momenta,
excluding a silent μ-dependent detailed
balance factor error.

## Reproduction and next physical priority

[Code](vector_neutrino_finite_bath_time.py),
[4 targeted tests](../tests/test_vector_neutrino_finite_bath_time.py),
[compact results](results/vector-neutrino-finite-bath-time-2026-11-22.json).

    python -m pytest -q tests/test_vector_neutrino_finite_bath_time.py
    python -m kinetics.vector_neutrino_finite_bath_time --tau 1000 --steps 100 --modes 96

**Remaining physics blocker:** the finite ν FD projection
presumes instantaneous rapid elastic redistribution
which is NOT derived from actual SM weak collisions.
It preserves combined number/energy and permits
an independent equilibrium comparison, but
cannot replace a full nonthermal ν momentum-source,
electron–positron/γ plasma, neutrino flavor
decoupling, Z′↔χ2χbar2/χ1χbar1↔2Z′,
physical gauge mediator couplings/thermal masses
and expanding FRW. The closed-model invariant
Nν+2NV is violated by additional SM ν number-
changing channels and must be dropped or replaced
with correct charges once those are included.

This is **not** a validated physical χ1–SM Γel,
not χ1↔χ2+φ Bose 86-eV relic evolution,
and yields no ΩDM, macroscopic BEC, frozen fφ(q),
CLASS or observational prediction.
