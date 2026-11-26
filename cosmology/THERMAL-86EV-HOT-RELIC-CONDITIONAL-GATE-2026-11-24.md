# CONDITIONAL hot relic / small-scale-structure screen: 86-eV real Bose scalar (2026-11-24)

## Purpose and scientific boundary

The open Dark-Matter research programme requires physical dark-sector Yukawa kinetics, finite Standard Model heat exchange, Friedmann expansion, frozen absolute fφ(q), Ωφ, CLASS perturbations and actual observations. The 2026-11-24 shared-event fully coupled χ1↔χ2+φ Bernstein-P2 *collision-only* prototype is not that cosmological solution. The project has **already verified this specific thermal benchmark with the actual CLASS solver** on 2026-11-21 (see [THERMAL-CLASS-VERIFIED.md](THERMAL-CLASS-VERIFIED.md) and [successful run #44890003670](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44890003670)). The NEW contribution of this stage is an independently integrated analytical RD+MD **mean path length** and a comparison to those existing CLASS results. We therefore test the same sharply defined conditional hypothesis rather than claiming a first thermal CLASS calculation:

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
| 86eV real BE, **same g*S,dec=55 as already verified CLASS** (our **rounded** nγ=410.7/cm³ gives Ωh²=0.11912165; prior CLASS analytical normalization Ωh²=0.11912939) | **0.01039704** | **6.62721 Mpc** |
| 5.7keV two-state FD all-DM, temperature adjusted to Ωc h²=0.120 | **3.96079×10−5** | **0.055475 Mpc** |

**Independent analytic-to-CLASS cross-check:** for the previously evaluated g*S,dec=55 thermal scalar, the new direct physical momentum integral gives mean free-streaming path **6.62721131 Mpc**. The pre-existing verified CLASS solver returned a linear matter-power amplitude-half crossing at z=5, k_half=**0.43275181165 Mpc⁻¹**, with <0.04% PSD-grid shift ([original independent CLASS report](THERMAL-CLASS-VERIFIED.md)); their product k_half × mean_path ≈**2.86794**. That product is **not a universal conversion between free-streaming and power-half scales**, but the two distinct computations consistently indicate large-scale thermal suppression. Existing [toy colder-PSD CLASS cases](COLDER-CLASS-RESULTS.md) and [hot-tail CLASS cases](VERIFIED-HOT-TAIL-CLASS.md) remain phenomenology, not Yukawa production outputs.

The 86-eV all-DM thermal scalar has a characteristic pre-equality streaming path **~120 times longer** than the 5.7-keV all-DM thermal fermion comparator. This is a physically meaningful **structure-formation warning**. [Iršič et al., Phys. Rev. D 109, 043511 (2024)](https://doi.org/10.1103/PhysRevD.109.043511) report a 5.7-keV (95% confidence) lower mass for their *specific thermal-relic WDM* Lyman-α analysis, relaxing to 4.1 keV for restricted scales/priors. **It is NOT statistically valid to call our 86-eV Bose scalar excluded at 95% simply by copying the fermion mass limit**; differing FD/BE momentum distributions, fractions, thermal histories, perturbations, and nuisance parameters require a **matched Bose CLASS/CAMB transfer function and observation-matched likelihood**. The required thermal transfer function already exists in the repo; only the observation likelihood and physically generated nonthermal Yukawa transfer remain uncalculated.

Nevertheless, a simple near-86-eV **fully thermal all-DM free-streaming hypothesis** is a poor cold-structure candidate: Ωh² can match while momenta remain too large to be treated as CDM without an actual transfer-function check. This is a conditional negative *screen*, not an exclusion of all 86-eV bosonic DM (nonthermal very cold spectra or genuine condensates are distinct hypotheses).

## Reproducibility, errors, and next action

Source: [independent analytic/numerical thermal screen](thermal_86ev_boson_screen.py).
Focused tests: [five checks](../tests/test_thermal_86ev_boson_screen.py).
Commands:

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_thermal_86ev_boson_screen.py
    OPENBLAS_NUM_THREADS=1 python -m cosmology.thermal_86ev_boson_screen

The momentum integration is validated against exact BE <x>=2.70117803 and FD <x>=3.15137437 and converged in order 48/96/128. Five local tests passed in 0.77 seconds; GitHub Scientific CI must be independently checked for new commits.

**Limitations:** g*S(1GeV)=61.75 is a rough step approximation; the actual freeze-out temperature and g*S(T) must follow the microscopic network and SM-plasma thermodynamics, including any dark↔SM entropy transfer. This free-streaming horizon is an unperturbed RD+MD diagnostic, not a power-spectrum transfer function. Assumes free decoupling from nearly zero scale factor (finite 1-GeV decoupling start gives negligible horizon correction) and one thermal non-condensed scalar with no later reheating. It is not an actual physical Yukawa freeze-out calculation, no CLASS-ready model, no observed-fit conclusion.

**NEXT_ACTION:** retain the **already verified** CLASS 86-eV thermal transfer as a conditional independent control, not a new discovery; first derive the **actual** nonthermal late-time Yukawa+SM+FRW fφ(q) and physical Ω without an inexhaustible reservoir. Pass THAT frozen PSD through the already functioning CLASS pipeline and then compare with the actual Lyman-α/CMB/galactic likelihood, without transplanting a fermion WDM mass bound.

**Input-normalization correction (CI negative result):** The independent free-streaming screening code deliberately uses the rounded PDG photon density nγ=410.7 cm⁻³. For g*Sdec=55 its own analytically computed Ωφ h² is **0.11912164989760488**, whereas the existing CLASS benchmark's more precise photon-density convention gives **0.11912939486124709** (relative difference **6.5013×10⁻⁵**). The first new cross-check test incorrectly assumed bitwise equality and failed Python 3.11/3.13 Scientific CI; it is corrected to verify each value against **its own declared inputs** and independently bound the relative normalization discrepancy. The mean free-streaming path uses the temperature ratio and is unaffected. This is an input precision difference, not a physical relic tension or an observation constraint.
