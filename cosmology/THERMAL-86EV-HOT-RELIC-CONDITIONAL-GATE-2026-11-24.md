# CONDITIONAL hot relic / small-scale-structure screen: 86-eV real Bose scalar (2026-11-24)

## Purpose and scientific boundary

The open Dark-Matter research programme requires physical dark-sector Yukawa kinetics, finite Standard Model heat exchange, Friedmann expansion, frozen absolute fφ(q), Ωφ, CLASS perturbations and actual observations. The 2026-11-24 shared-event fully coupled χ1↔χ2+φ Bernstein-P2 *collision-only* prototype is not that cosmological solution. Before expensive FRW+CLASS, it is scientifically useful to screen an analytically solvable **alternate, sharply defined hypothesis**:

> What if the 86-eV real scalar had a **zero-chemical-potential thermal Bose-Einstein spectrum**, decoupled while relativistic from an SM bath around 1 GeV, and subsequently simply free-streamed with conserved comoving entropy, no condensate and no extra dark reheating?

This is **not** an inferred property of the existing vacuum-Yukawa model. It is an independent conditional benchmark / possible falsification gate. The numerical values are reproducible with the script and use explicit observable cosmological inputs, not a hidden reservoir in the **evolving** Yukawa solver.

Observational context: [Planck Collaboration, 2018 VI](https://doi.org/10.1051/0004-6361/201833910) quotes Ωc h²≈0.120 and Ωm≈0.315, H0≈67.4 km/s/Mpc. [PDG astrophysical constants](https://pdg.lbl.gov/2024/download/db2024.pdf) provide Tcmb=2.7255K, nγ≈410.7 cm−3, ρcrit/h²≈1.053672×10⁴ eV cm−3, and zeq≈3402. A simple step-function SM thermal-degrees count for a ~1-GeV benchmark is g*S,dec≈61.75; this is NOT a measured decoupling temperature or a replacement for smooth lattice-QCD g*S(T).

## Independent relativistic-decoupling abundance

For one real bosonic degree with zero μ, no BEC and a conserved BE spectrum after relativistic decoupling,

    (Tφ/Tγ)^3 = g*S,SM(today)/g*S,SM(dec).
    nφ0 = (1/2) nγ0 (Tφ/Tγ)^3.
    Ωφ h² = mφ nφ0 / (ρcrit/h²).

Take g*S,SM(today)=3.909, g*S,SM(dec)=61.75, mφ=86eV. This yields:

| Conditional thermal observable | Analytic result |
|---|---:|
| Tφ0/Tγ0 | **0.3985439641** |
| Tφ0 | 9.360419474e−5 eV |
| nφ0 | **12.99940324 cm−3** |
| Ωφ h² | **0.106100255** |
| ΔNeff (when φ is still relativistic and decoupled) | **0.05554485** |
| g*S,dec required to match Ωc h²=0.120 at 86eV | **54.59742287** |

The abundance can fall near the observed dark-matter density by choosing a plausible but **currently unproved** entropy/decoupling history. An abundance match alone does NOT establish a viable cold dark-matter model. The true χ+SM+φ kinetics may be nonthermal, may violate this hot-freezeout ansatz, or may involve an independent zero-mode.

## Direct integral over the entire Bose momentum distribution

The Bose number-weighted comoving momentum probability is

    P_B(x)dx = [x²/(2ζ(3)(exp(x)−1))] dx, x=q/Tφ0.
    <x>_B = π⁴/(30ζ(3)) = 2.7011780329.

For a collisionless particle with fixed comoving momentum q, physical velocity at scale factor a is v(q,a)/c=q/sqrt(q²+(mφ a)²). In a background with matter Ωm=.315 and equality zeq=3402, set Ωr=Ωm/(1+zeq), aeq=(1+zeq)^−1. The **mean traveled comoving path** before equality in a free-streaming RD+MD approximation is

    <λ_fs(aeq)> = (c/H0) ∫_0^aeq da
                   [∫ P_B(q) q/sqrt(q²+(mφ a)²) dq]
                   /sqrt(Ωr+Ωm a).

This is computed independently from the Yukawa collision solver by exact momentum distribution integration using positive Gauss-Laguerre quadrature (orders 48/96/128) and adaptive SciPy quad_vec on scale factor. It is a **characteristic free-streaming path**, **not** the linear Boltzmann/CLASS transfer function, damping half-mode length, or a direct Lyman-α likelihood.

| Thermal benchmark | <v/c> at equality | <λ_fs> to equality |
|---|---:|---:|
| 86eV real BE, g*S,dec=61.75 (Ωh²=0.1061) | **0.01000359** | **6.42844 Mpc** |
| 86eV real BE, temperature adjusted to Ωc h²=0.120 | **0.01042252** | **6.64002 Mpc** |
| 5.7keV two-state FD all-DM, temperature adjusted to Ωc h²=0.120 | **3.96079×10−5** | **0.055475 Mpc** |

The 86-eV all-DM thermal scalar has a characteristic pre-equality streaming path **~120 times longer** than the 5.7-keV all-DM thermal fermion comparator. This is a physically meaningful **structure-formation warning**. [Iršič et al., Phys. Rev. D 109, 043511 (2024)](https://doi.org/10.1103/PhysRevD.109.043511) report a 5.7-keV (95% confidence) lower mass for their *specific thermal-relic WDM* Lyman-α analysis, relaxing to 4.1 keV for restricted scales/priors. **It is NOT statistically valid to call our 86-eV Bose scalar excluded at 95% simply by copying the fermion mass limit**; differing FD/BE momentum distributions, fractions, thermal histories, perturbations, and nuisance parameters require a Bose CLASS/CAMB transfer function and observation-matched likelihood.

Nevertheless, a simple near-86-eV **fully thermal all-DM free-streaming hypothesis** is a poor cold-structure candidate: Ωh² can match while momenta remain too large to be treated as CDM without an actual transfer-function check. This is a conditional negative *screen*, not an exclusion of all 86-eV bosonic DM (nonthermal very cold spectra or genuine condensates are distinct hypotheses).

## Reproducibility, errors, and next action

Source: [independent analytic/numerical thermal screen](thermal_86ev_boson_screen.py).
Focused tests: [five checks](../tests/test_thermal_86ev_boson_screen.py).
Commands:

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_thermal_86ev_boson_screen.py
    OPENBLAS_NUM_THREADS=1 python -m cosmology.thermal_86ev_boson_screen

The momentum integration is validated against exact BE <x>=2.70117803 and FD <x>=3.15137437 and converged in order 48/96/128. Five local tests passed in 0.77 seconds; GitHub Scientific CI must be independently checked for new commits.

**Limitations:** g*S(1GeV)=61.75 is a rough step approximation; the actual freeze-out temperature and g*S(T) must follow the microscopic network and SM-plasma thermodynamics, including any dark↔SM entropy transfer. This free-streaming horizon is an unperturbed RD+MD diagnostic, not a power-spectrum transfer function. Assumes free decoupling from nearly zero scale factor (finite 1-GeV decoupling start gives negligible horizon correction) and one thermal non-condensed scalar with no later reheating. It is not an actual physical Yukawa freeze-out calculation, no CLASS-ready model, no observed-fit conclusion.

**NEXT_ACTION:** preserve this conditional thermal relic as a hard *screening benchmark* for any future coupled SM+FRW fφ(q); test whether physical collision-time evolution creates a much colder spectrum (or a consistently derived zero mode), calculate normalised Ω and ΔNeff, then generate CLASS transfer functions and compare with the published actual Lyman-α likelihood / galaxy small-scale data, without transplanting a fermion WDM bound.
