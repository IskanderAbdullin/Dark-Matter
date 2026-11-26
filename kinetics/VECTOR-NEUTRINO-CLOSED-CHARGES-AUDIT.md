# Conserved ν-pair number prevents zero-chemical-potential Z′ equilibrium in a closed ν+Z′ quantum model (2026-11-22)

## The scientific question and precise scope

Previous newly derived microscopic moving Z′↔ννbar finite-T Pauli/Bose source shows **fast initial inverse production** relative to an assumed H, under an arbitrary contact G and perturbative dark coupling. But rapid 1↔2 collision *by itself* does not imply the mediator reaches a canonical zero-chemical-potential BE distribution with no energy or particle-number backreaction. What is its actual **closed-system equilibrium** if ONLY the two-active-flavor reactions Z′↔νeνbar_e and Z′↔νμνbar_μ operate, and initially Z′ is absent?

This audit solves the constrained maximum-entropy **equilibrium**, not time evolution, in an exactly specified reduced system: two active massless neutrino flavors + antineutrinos (4 individual single-helicity states), one broken-phase massive spin-one Z′ (3 polarizations), fixed physical volume, common eventual T, no electron/positron/photon/χ couplings, no expansion, ν–ν elastic redistributions assumed sufficient for kinetic equilibration, and zero lepton-flavor **asymmetry**. Both ν and antineutrino have a common effective *particle-number* chemical potential μν, **not** a lepton asymmetry; the boson has μV=2μν from reaction equilibrium. This closure is only a contrast to imposing μν=μV=0. The real SM violates the reduced number invariant through weak pair reactions, so the result is NOT a physical standard-cosmology abundance.

Pre-registered checks: exact total energy and genuine collision number charge relative residual <1e−11, entropy must increase, radial n128→256 T and μ/T agreement within 1e−8; independent zero-mass-ratio fixed-three-polarization fugacity series relative T within 4e−5 and μ/T within 2e−4 at mV/T0=.03. Do not declare condensate or physical Neff.

## Exact conserved quantity and chemical equilibrium

Each elementary reaction removes one boson and creates exactly ν_f+νbar_f, or its inverse. Therefore when **all other number-changing channels are absent**,

    N_total = N_ν+N_νbar+2 N_V

is conserved, as is total physical energy E_total=E_ν+E_νbar+E_V. Neutrino lepton-flavor **net** charges Nν−Nνbar separately remain zero; individual neutrino and vector numbers are NOT conserved. This is independent of electric-charge conservation in the e± bath and must not be conflated with it.

For each of ν and antiν, the FD distribution is
fν(q)=[exp(q−y)+1]^−1 with y=μν/T.
The mediator BE distribution is
bV(q)=[exp(sqrt(q²+x²)−2y)−1]^−1, x=mV/T,
and μV=2μν. For the reference branch with negative μν there is no Bose chemical-potential singularity.

With dimensionless q=p/T, define

    nν=4T³/(2π²) ∫dq q² fν(q),
    ρν=4T⁴/(2π²) ∫dq q³ fν(q),
    nV=3T³/(2π²) ∫dq q² bV(q),
    ρV=3T⁴/(2π²) ∫dq q² sqrt(q²+x²)bV(q).

The unknown T and y solve

    nν(T,y)+2nV(T,2y)=nν(T0,0),
    ρν(T,y)+ρV(T,2y)=ρν(T0,0),

for initial thermal neutrinos at T0, μν=0 and **exactly empty** mediator bV=0.

The total equilibrium entropy is calculated independently:
Sν=[(4/3)ρν−μν nν]/T,
SV=3/(2π²)∫dp p²[(1+bV)ln(1+bV)−bV ln bV].
A physical equilibrium must increase entropy from the initial vector-empty state.

## Numerical solution and negative counterexample to naive thermal reset

Use energy mass unit T0=1, mV/T0=.03,.1,.3;
Gaussian p/T quadrature n=256 on [0,60],
two-dimensional root solve for (ln(T/T0), μν/T), and
**two independent conserved residuals**.

| mV/T0 | T_final/T0 | μν/Tfinal | nV_final / nV_BE(μ=0,T0) | nν_final/nν_initial | S_final/S_initial |
|---:|---:|---:|---:|---:|---:|
| .03 | 1.24600008 | −1.13589978 | **0.16828877** | 0.66371856 | 1.07289863 |
| .10 | 1.24565840 | −1.13489124 | **0.16932222** | 0.66381836 | 1.07283951 |
| .30 | 1.24276282 | −1.12629010 | **0.17555559** | 0.66469279 | 1.07233249 |

For mV/T0=.03 the final vector number equals
0.16814072×the initial ν+νbar number; hence
nν_final/nν0+2 nV_final/nν0
=0.66371856+2×0.16814072=1
to double precision. The total energy residual is
≤5e−16 relative and entropy increases ~7.29%.
n128→256 shifts T/T0 and μ/T at ≲1e−11 relative
and the quoted vector ratio by <2e−7 relative
for mV/T0=.03; the latter is more IR-sensitive.

**Explicit negative result:** even if the microscopic
1↔2 process became arbitrarily fast, its closed-system
chemical equilibrium would NOT simply be μν=μV=0.
It instead develops a *negative common particle-number
chemical potential*, with a vector number only ~17%
of the μV=0 BE reference at the original T0, due to
the conserved Nν+2NV and energy. A model that
instantaneously resets the neutrino and mediator
distributions to independent zero-μ thermal states
would **create particles/energy from outside the
declared closed system**, regardless of green CI.

The temperature rises from T0 to ~1.246T0 in this
ideal constrained equilibrium even though some energy
moves to the mediator; the negative ν chemical
potential reduces neutrino number/energy at fixed T,
and conservation fixes the final coupled T.
This is not photon reheating, a cosmic neutrino
history, or a physical claim of Z′ thermalization.

## Independent nonquadrature reference

For a **formal mV/T→0 limit retaining three massive-vector
Proca polarizations** (NOT the restored-gauge physical
two-polarization phase), the moments have independent
absolutely convergent fugacity series at y<0:

    nν=4T³/π² sum_{k>=1} (−1)^(k+1)e^(ky)/k³,
    ρν=12T⁴/π² sum_{k>=1} (−1)^(k+1)e^(ky)/k⁴,
    nV=3T³/π² sum_{k>=1} e^(2ky)/k³,
    ρV=9T⁴/π² sum_{k>=1} e^(2ky)/k⁴.

At initial y=0, exact nν0=3ζ(3)/π², ρν0=7π²/60.
Solving this **independently of the massive p/T
quadrature** gives

    T_final/T0 = 1.24603411773957,
    y = −1.13600012538095,
    nV/nν_initial = 0.16814566001216.

The mV/T0=.03 finite-mass result approaches this
formal limit at the predicted <4e−5 fractional
temperature and <2e−4 chemical-potential tolerance.
Neither this formal zero-mass continuation nor a
massive-vector 3-polarization equilibrium can be used
across unknown U(1)′ symmetry restoration without
including its gauge-Higgs sector.

## Reproduction and next physical step

[Executable code](vector_neutrino_equilibrium_charges.py),
[4 direct/independent tests](../tests/test_vector_neutrino_equilibrium_charges.py),
[compact results](results/vector-neutrino-equilibrium-charges-2026-11-22.json).

    python -m pytest -q tests/test_vector_neutrino_equilibrium_charges.py
    python -m kinetics.vector_neutrino_equilibrium_charges

Next priority: use the already derived *exact*
mode-resolved Pauli/Bose Z′↔νν collision operator
to integrate finite bV(p,t) **and** finite ν
temperature/chemical potential evolving under both
Nν+2NV and total energy conservation, with positive
entropy and independently converged dt/q/cutoff;
then turn on actual SM weak number-changing
ννbar↔e+e−, Z′↔χ_iχbar_i and χ1χbar1↔2Z′,
FRW/H and proper gauge-sector thermal masses.

**Scope:** no full ν phase-space strong collision or
actual weak elastic thermalization times, no cosmological
Z′ abundance, finite SM e±γ, χ1/χ2 nonthermal kinetics
or χ1↔χ2+φ Bose 86-eV freezeout.
This is NOT a DM exclusion, macroscopic BEC,
absolute Ωφ or CLASS-ready prediction.
