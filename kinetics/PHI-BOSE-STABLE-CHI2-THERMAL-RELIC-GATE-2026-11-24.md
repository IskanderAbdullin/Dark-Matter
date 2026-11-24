# Stable Dirac chi2 partner: analytic p-wave Yukawa and exact thermal symmetric relic negative (2026-11-24 UTC)

## Why this is a necessary gate for the global Bose dark matter goal

The real-phi renormalizable Yukawa y phi (bar chi2 chi1 + bar chi1 chi2), with m1=100 MeV, m2=20 MeV, mphi=86 eV, preserves the **same net Dirac U1** carried by chi1 and chi2. A consistent candidate vector mediator has q1=q2 under this U1; kinetic mixing gives electrons a vector coupling but does NOT erase dark fermion number. In this explicitly minimal field content, chi1 can decay to chi2+phi, but the lightest charged fermion **chi2 and its antiparticle cannot decay into neutral phi alone**. Therefore total dark density must account for **Omega_phi + Omega_chi2+anti2**, not only Omega_phi.

This stage asks a falsifiable, restricted question: IF a symmetric thermal chi2/antichi2 bath is ever populated and remains in MB thermal contact with a thermal phi bath until chemical freezeout, and IF the only appreciable number-changing late channels are chi2+antichi2 -> phi+phi via virtual chi1 plus the offshell vector into e+e−, does y=0.001 allow the stable 20-MeV component to be subdominant? **No, very strongly, within those assumptions.** This is not an exclusion of all U1 bosonic dark matter histories; a never-thermalized sector, new interactions, chemical asymmetry, nonthermal Bose factors, or different masses/couplings require new modeling.

## Analytic p-wave derived from the full physical t/u tree amplitude

For two **identical real final scalars** (include final 1/2!), the already independently implemented spin-average is the full Dirac t+u exchange of virtual chi1 of mass M. Expanding the exact 4×4 Dirac gamma trace for chi2 mass m, negligible mphi/m≈4.3e−6, and integrating CM angle gives the exact coefficient in the v_Moller=2p/E convention:

    sigma(chi2+antichi2 -> phi+phi)*v_Moller
      = b * v_Moller^2 + O(v_Moller^4, mphi^2/m^2)

    b = y^4*(3 M^6 + 6 M^5 m + 9 M^4 m^2
             + 8 M^3 m^3 + 5 M^2 m^4
             + 2 M m^5 + 3 m^6)
             /[96 pi (M^2+m^2)^4].

For M=.1GeV, m=.02GeV, y=.001, the coefficient is **b=1.3130827280563e−12 GeV^-2**. There is NO s-wave term: at zero relative momentum the fully summed t/u spin-averaged amplitude vanishes. The analytical expression is validated by a genuinely independent complex gamma-matrix trace of the two physical propagators with 64 Gaussian CM angles and by the older vectorized invariant-trace sigma(s) at beta=.02,.05,.1,.2. The remaining discrepancy shrinks as beta².

In the same portal candidate mV=10GeV,gD=.1,epsilon=1e−5, direct offshell chi2+antichi2 -> e+e− has s-wave

    a_e = (gD epsilon e)^2*m^2 /
          [pi ((4m^2-mV^2)^2 + mV^2 GammaV^2)]

which equals **1.167613445795e−21 GeV^-2**, utterly negligible relative to the Yukawa p-wave 6b/x around x=m/T≈7–15 for this specific benchmark. A REAL on-shell mediator resonance can dominate at T~1GeV but is not kinematically accessed by nonrelativistic chi2 at T~few MeV; do not use its GeV forward/H number as the late stable partner annihilation rate.

## Conditional radiation-dominated abundance calculation

Define Y=n_chi2/s for ONE Dirac charge sector (spin g=2), and symmetric n_antichi2=n_chi2. Impose **constant** g*energy=g*entropy=10.75, external RD H(T), MB n_eq=g m² T K2(m/T)/(2pi²) and Y(x=3)=Y_eq(3), with x=m/T. Integrate the controlled one-sector Boltzmann equation

    dY/dx = -s/(H x) <sigma v> (Y^2-Yeq^2),
    Omega_(chi2+antichi2) h² = 2*m*(s0/[rho_crit/h²])*Y_infinity,

with s0=2891.2cm^-3 and rho_crit/h²=1.05375e−5 GeV/cm³. No hidden normalization to the observed DM abundance and no assumption that phi makes all DM. This is an explicitly **simplified thermal branch**, not the full coupled chi1/chi2/phi/SM Friedmann system; phi is assumed to have negligible chemical potential and to provide equilibrium inverse reactions.

There are TWO distinct numerical estimators. An asymptotic NR p-wave <sigma v>=a_e+6b/x throughout x≥3 yields **Omega h²≈761.289**, but this is NOT fully controlled at the actual x~7–10 chemical decoupling. We therefore independently thermally average the **full** vacuum χ2 pair t/u cross section with the relativistic Gondolo–Gelmini Bessel integral, shifting u=sqrt(s)/T to u=2x+t to cancel dangerous exp(-2x) factors analytically. Its ratio to naive p-wave is:

| x=m/T | True relativistic MB <sigma v> / (6b/x) |
|---:|---:|
| 3 | 0.481333 |
| 8 | 0.748759 |
| 20 | 0.888606 |
| 100 | 0.976404 |

Thus the naive p-wave approximation **overestimates** actual thermal annihilation strength by ~34% at x8 (when compared to the actual exact rate) and understates the residual relic. Integrating with the exact thermal GG rate from x3 to x100, then an explicitly declared p-wave/asymptotic empty-equilibrium tail from x100 to infinity gives:

| logarithmic exact thermal interpolation samples | Omega_(chi2+antichi2) h² |
|---:|---:|
| 16 | 879.95226363 |
| 32 | 879.95844334 |
| 64 | **879.95876900** |

The last two differ by 3.7e−7 relative; this measures interpolation resolution for the stated fixed-g* ODE, NOT cosmological model error. At y=.001 the outcome is approx **7333 times** the commonly quoted Planck base-LCDM cold DM density **Omega_c h²≈0.120±0.001** ([Planck Collaboration 2020, A&A 641 A6](https://doi.org/10.1051/0004-6361/201833910)). The Planck number is a target for a physically correct late CDM sector, **not** an observed test of this restricted model alone.

An illustrative inversion of THIS toy ODE finds that y would have to rise to about **0.0135** merely to make the stable chi2+antichi2 component alone match 0.120. If phi is also to be dominant, chi2 must be more suppressed still. But y changes BOTH the χ1↔χ2φ decay and χχbar↔φφ collision rates and their full thermal/FRW history; 0.0135 is **not a fitted, experimentally allowed, or complete physical solution**. Even the exact pair thermal test assumes a vacuum mediator, no medium/thermal masses and simple MB bath.

## Key missing physics and logically permitted conclusion

- The real 5-population FD/BE source has four channels and a finite SM heat exchange; this reduced equation omits chi1 replenishment, nonthermal excitation, Pauli/Bose factors at freezeout, evolving phi temperature and chemical potential, condensate and null modes, actual g*energy/entropy changes across e± reheating, thermal mass, width, realistic SM species and Hubble feedback. Thus it is not a strict exclusion of the full original candidate.
- However, **a stable 20-MeV fermion partner cannot simply be ignored while claiming a full 86-eV bosonic dark matter model**. Its abundance and phase-space distribution must be independently computed, or a validated interaction that depletes/decays it must be added and retested under the same gauge/Z2 constraints.
- The suggested vector A' has other open SM charged decays and hadronic uncertainties; the electron-only annihilation and external R(s) gate are recorded in the companion vector report. When a dynamic A' is included explicitly, avoid double counting real-intermediate-state resonance and full 2↔2 diagrams.

**Reproduction:** [analytic/ODE/GG source](phi_bose_stable_chi2_relic_gate.py), [independent Dirac trace, Bessel and time/thermal mesh tests](../tests/test_phi_bose_stable_chi2_relic_gate.py), [permanent machine data](results/phi-bose-stable-chi2-thermal-gate-2026-11-24.json), previous independent [exact pair sigma](phi_bose_pair_independent_vectorized.py). Five-field focused CI #45452412200 success for full Gamma/thermal-grid validation; later exact-y-target test run to be checked separately.

**NEXT ACTION:** the scientifically highest-value physical gate is now to evolve BOTH stable chi2 and phi through an actual coupled, quantum FD/BE dark+SM finite heat/FRW system including reverse pair and number-changing SM processes; compare absolute Omega_chi2+anti2+phi vs measured OmegaDM with true g*(T) and uncertainty. If chi2 overcloses even across credible physical variations, identify an explicitly excluded parameter subregion; do not convert this toy result into universal DM exclusion.
