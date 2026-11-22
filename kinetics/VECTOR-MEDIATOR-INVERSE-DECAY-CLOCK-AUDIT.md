# Conditional vector-neutrino decay and Hubble gate: fixed contact product is not a virtual-only plasma (2026-11-22)

## Scope, physical prerequisites and pre-declared numerical controls

The candidate vectorlike χ1/χ2 gauged U(1)_(L_e−L_μ) model with a separate massive Z′ has a **necessary tree-level decay channel** Z′→νe νbar_e,νμ νbar_μ (two massless active flavors) even where Z′→e+e− is kinematically closed. Prior thermal χ+e/ν scattering benchmarks used arbitrary zero-momentum coupling G=0.001 m_e^-2 and mediator mV=.03,.1,.3 m_e, but never computed *real Z′ population or decay*.

A single contact coefficient G=(gχ qχ)g′/mV² fixes only the **product** of the dark and lepton couplings. To obtain a conditional lower estimate of the required lepton coupling, impose an explicit *conventional perturbative bound* |gχqχ|≤√(4π) on the effective dark vector coupling, assume same vector mass and contact G, and choose the MAXIMUM allowed dark coupling. This minimizes g′ and its neutrino width within this chosen convention. It does **not** enforce experimental limits, gauge symmetry restoration, scalar masses, running coupling or even compatibility with the eventual χ1,χ2,φ production history.

Validation: independently evaluate the Maxwell–Boltzmann time-dilated decay factor via Bessel K1/K2 and a separate momentum integral, relative difference <1e−9 for x=mV/T=.03,.1,.3,1,3. Compare decay clock with a **stated reference** radiation H using gstar=10.75 and non-reduced Mpl=1.2209e19 GeV; never identify a ratio with an observed thermalization rate without the quantum inverse source.

## Tree neutrino width, time dilation and Pauli bound

With two massless active ν flavor-antiflavor channels and left U(1)′ coupling g′,

    Gamma_νν,rest = g′² mV/(12π).

For a moving unpolarized vector E=sqrt(p²+mV²), vacuum decay clock is Gamma_rest mV/E. Averaging this over an illustrative **Maxwell–Boltzmann distribution of existing vectors**, but **without** neutrino Pauli final-state blocking and inverse scattering:

    <Gamma>_MB = Gamma_rest K1(x)/K2(x),
    x=mV/T,

where K1/K2 is independently equal to

    [integral dp p² exp(-E/T) mV/E] /
    [integral dp p² exp(-E/T)].

This is NOT the thermal quantum chemical relaxation eigenvalue. If the thermal ν and antiν distributions are FD with zero chemical potential and positive energies, each f_FD≤1/2, so each final Pauli factor (1−f)≥1/2. Thus, **conditional on the same MB momentum distribution of existing vectors and no other plasma effects**,

    <Gamma_νν,loss>_MB ≥ <Gamma_νν,vacuum>_MB/4.

This is a rigorous lower bound for **loss of an existing vector via νν**, not for *net equilibration*, because inverse νν→Z′ is absent. If vectors are initially absent, one must integrate inverse production instead of applying a decay lifetime to an unpopulated species.

## Quantitative thermal-vs-expansion diagnostic for the exploratory G

At T=m_e=0.51099895 MeV, with radiation gstar_rho artificially fixed to 10.75 (no extra Z′/dark backreaction), H≈1.164051626×10⁻²⁵ GeV. Choose |gχqχ|=√(4π)=3.5449077, the boundary of the stated perturbative convention.

| mV/me | g′ required by G mV² | K1/K2 | Γνν,MB/H | Existing-MB-vector FD Pauli-loss lower /H |
|---:|---:|---:|---:|---:|
| .03 | 2.5388531e−7 | 0.01497554 | 3.3721e3 | ≥843 |
| .10 | 2.8209479e−6 | 0.04939172 | 4.5768e6 | ≥1.14e6 |
| .30 | 2.5388531e−5 | 0.14053290 | 3.1644e9 | ≥7.91e8 |

The ratio scales as g′², so smaller perturbative |gχqχ| for the **same fixed G** would require larger g′ and increase these vacuum/MB clocks. This is a **specific, conditional strong inconsistency warning** for interpreting the arbitrary earlier contact benchmarks as a weakly coupled virtual-only neutrino plasma: the gauge mediator has an open very fast neutrino loss/production channel relative to expansion at those masses and couplings. Nonetheless, the ratio alone cannot certify a thermal vector population or observational exclusion, because exact fermion blocking, stimulated vector inverse production, ν chemical potentials, actual Tν(T), and mV(T) have not been computed.

**Negative distinction:** Γνν,MB/H and its Pauli lower bound are not Γ_el/H and not proof of the 86-eV Bose-production mechanism. One may reduce G or invoke a different mediator, scalar field content, higher-temperature symmetry restoration, nonthermal ν bath or other UV completion—each is a DIFFERENT physical parameter/model and must be audited, not silently substituted into the earlier 5D χ1-only collision.

## Independent source, reproduction and next fundamental gate

[Conditional code](vector_mediator_inverse_decay_gate.py) ·
[8 regression checks](../tests/test_vector_mediator_inverse_decay_gate.py) ·
[compact results](results/vector-mediator-inverse-decay-2026-11-22.json).

    python -m pytest -q tests/test_vector_mediator_inverse_decay_gate.py
    python -m kinetics.vector_mediator_inverse_decay_gate

Literature scope: [N_eff constraints on light U(1) gauge bosons, Phys Rev D 110, 075032 (2024)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.110.075032), [light vector effect on neutrino decoupling and energy transfer, EPJC 2024](https://link.springer.com/article/10.1140/epjc/s10052-024-13220-8). They motivate a quantitative BBN/CMB calculation but DO NOT validate the arbitrary G/mV or establish an exclusion here.

**Immediate next P0:** implement full finite-temperature **mode-resolved Z′↔νν** collision, with real moving-vector kinematics, neutrino FD blocking, inverse νν annihilation and Z′ Bose enhancement, test thermal detailed balance/energy/number/entropy, independently integrate vector production and compare to MB estimates and H. Then add actual Z′→χ2χbar2 and χ1χbar1↔2Z′ rates, gauge-S sector and constraints, evolving χ1/χ2/real φ plus e±γν/Z′ to FRW and freezeout with no hidden bath. No physical χ Γ_el, 86-eV ΩDM, BEC or CLASS-ready prediction exists.
