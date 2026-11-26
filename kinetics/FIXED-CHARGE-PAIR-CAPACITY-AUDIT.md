# Fixed-charge thermal susceptibility of the massive e± bath (2026-11-22)

## Motivation and pre-declared validation gate

The previous ideal e±–γ EOS demonstrates finite heat and the
necessity of nonzero electron chemical potential at low T. However,
its simple C=∂ρ/∂T at **μe/T held fixed** is generally *not*
the heat capacity available in a finite-volume 2→2 χ–e collision,
which conserves net electric charge. Derive the fixed-charge capacity
using the exact FD susceptibility matrix, then independently check
by finite differences with **the same physical charge density** at
T±δT. Numerical tolerance for T/me=.05,.1,1,.025:
relative deviation ≤3×10⁻⁶ for h=10⁻⁴.

The initial background is an **illustrative** net charge density
nQ=ηnγ(Tinitial), η=6×10⁻¹⁰. At a fixed collision volume, nQ
is held constant when varying T; η is **not** held constant in
the derivative, since nγ∝T³. No baryon dynamics, cosmological
η(T), photon reheating, annihilation rates or ν/QED physics
is inferred.

## Derivation

Let x=me/T, q=p/T, ε=√(q²+x²), y=μe−/T=−μe+/T, and
f−=1/[exp(ε−y)+1], f+=1/[exp(ε+y)+1].
At a temperature T with a specified nQ=n_e−−n_e+,

    nQ/T³ = (1/π²) ∫dq q²(f−−f+),
    ρe±/T⁴ = (1/π²) ∫dq q² ε(f−+f+).

Set a=f−(1−f−), b=f+(1−f+). Three positive/signed
dimensionless thermodynamic susceptibilities are

    B=(1/π²)∫dq q² ε² (a+b),
    D=(1/π²)∫dq q² ε   (a−b),
    S=(1/π²)∫dq q²     (a+b).

Then, for a fixed volume,

    (∂ρe±/∂T)|y = T³ B,
    (∂ρe±/∂y)|T = T⁴ D,
    (∂nQ/∂T)|y = T² D,
    (∂nQ/∂y)|T = T³ S,

and the properly **charge-constrained** electron-pair heat
capacity is exactly

    C_e±|nQ /T³ = B − D²/S.

The subtracted term is a Schur complement, not an arbitrary
numerical correction. By Cauchy–Schwarz, 0≤D²≤BS; hence
C_e±|nQ≥0. The photon heat is added separately,
Cγ/T³=4π²/15. This is a local equilibrium susceptibility
for a *charged finite-volume* EM bath, not a rate or solution
of the quantum Boltzmann collision equation. The formula
assumes pair chemistry fast enough to maintain the two
FD functions with opposite y; this must be checked
microscopically and may fail after annihilation freeze-out.

## Numeric findings and negative result

p/T Gauss–Legendre n256 on [0,48], same massive FD
integrals as the independent EOS audit; solve neutrality
using stable sinh(y)/(cosh ε+cosh y) to avoid tiny-charge
cancellation. Electron/positron charge is constrained by
ηnγ at the *initial* temperature. For an independent
finite-difference check, solve y again at T±10⁻⁴T while
holding the original *physical density* nQ unchanged.

| T/me | y=μe−/T | B: C_e± at fixed y /T³ | C_e± at fixed charge /T³ | ratio C_Q/B |
|---|---:|---:|---:|---:|
| 1 | 5.08190e−10 | 4.428887902 | 4.428887902 | 1 |
| 0.1 | 3.35315e−7 | 0.06015415758 | 0.06015415758 | ≈1 |
| 0.05 | 0.002848678 | 2.399907705e−5 | 2.399888300e−5 | 0.999991914 |
| 0.025 | 13.83768468 | 2.524986517e−7 | 2.324267435e−10 | 0.000920507 |
| 0.02 | 23.51213620 | 3.884166249e−7 | 2.298649660e−10 | 0.000591800 |

**Negative thermodynamic result:** at T/me=.025, reusing
the fixed-y electron heat as if it were a fixed-charge heat
overestimates the electron capacity by about **1,086×**
in this illustrative charged gas. Crucially, photons still
dominate the *total* EM heat at such temperatures, so this
does not imply a corresponding factor in total SM heat.

Independent fixed-nQ finite-difference derivative compared with
B−D²/S: relative differences 1.04e−8 (T/me=1),
1.34e−7 (.1), 5.84e−7 (.05), and 7.55e−7 (.025).
At .025 the leading nonrelativistic residual-electron
classical estimate (3/2)nQ/T³=2.19229e−10 is only
~6.0% below the full finite-mass FD result, a further
asymptotic consistency check.

## Reproduction and limits

Executable source [fixed_charge_pair_capacity.py](fixed_charge_pair_capacity.py),
regression [test_fixed_charge_pair_capacity.py](../tests/test_fixed_charge_pair_capacity.py),
[compact machine-readable results](results/fixed-charge-pair-capacity-2026-11-22.json).

    python -m pytest -q tests/test_fixed_charge_pair_capacity.py
    python -m kinetics.fixed_charge_pair_capacity --T-over-me 0.025

No e+e−↔γγ collision rate, chemical freeze-out or photons/ν
energy transfer has yet been derived. Electrons actually
freeze into a constrained network with baryon charge and finite
chemical reaction times; one can neither choose the μ=0
equilibrium capacity at all T nor assign the fixed-charge
fast-chemistry capacity without proving reaction timescales.
The χe vector-contact EFT and e/ν gauge completion are
conditional; no full nonthermal 2→2+1↔2 evolution, physical
γel, BEC, Ωφ, frozen 86-eV f(q) or CLASS result follows.

**Next P0:** construct a charge/energy-conserving finite
e−/e+/γ/ν bath with explicit photon pair production/
annihilation and weak channel collision kernels (including
Bose/Pauli); compare full Boltzmann chemical relaxation
to both fast-chemistry and frozen-pair heat limits.
Then couple the validated χe nonlinear collision heat
without hidden thermal reset and bound EFT validity.
