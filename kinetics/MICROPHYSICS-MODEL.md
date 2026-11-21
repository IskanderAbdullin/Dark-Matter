# Microscopic portal pilot: explicit Yukawa decay and electron-scattering EFT

**Status:** conditional tree-level, finite-temperature collision operator **not yet derived**, not a validated cosmological 86-eV model. Code: [microphysics_portal.py](microphysics_portal.py); [unit and independent t-integral tests](../tests/test_microphysics_portal.py). Uses natural units with masses, T and rates in GeV. This is phase A2/A3 partial work in [research specification](CONTINUOUS-RESEARCH-SPEC.md).

## Falsifiable target and competing possibilities

We need replace the old arbitrary elastic contact rate `gamma_el=10^5` in [BEST+SM](BEST-SM-COUPLED-FRW.md) with a microscopic scattering operator, **without** assuming any particular 86-eV relic has been produced. First prove a concrete theory has the 1↔2 and SM contact vertices, their dimensions and thresholds; then test the full Fermi/Bose collision integrals, thermal rates, SM heat and freezeout. A high rest-frame scattering frequency alone is **not** evidence of thermalisation. If the portal's additional 2→2/number-changing processes invalidate the intended system, report the failure and change the model explicitly.

## Minimal low-energy benchmark (not a UV-complete SM model)

Fields: Dirac fermions chi1 and chi2, real scalar phi; exact Z2 with chi1 and phi odd, chi2 and SM even. For `m_phi < m_chi1`, phi is lightest Z2-odd particle and stable *within this symmetry*. A conserved Dirac charge can be assigned to chi1 and chi2. The interaction is

```text
L_int = -y phi (bar(chi2) chi1 + bar(chi1) chi2)
        + G (bar(chi1) gamma_mu chi1)(bar(e) gamma^mu e).
```

`y` dimensionless, `G=C/Lambda^2` with energy dimension -2. The electron vector contact describes the **broken-electroweak** low-energy component; at temperatures above EWSB it must be embedded in gauge-invariant operators such as `(bar chi1 gamma chi1)(bar L_e gamma L_e + bar e_R gamma e_R)`, necessarily carrying neutrino channels. A heavy mediator UV completion must also be supplied if the typical momentum transfer approaches Lambda. Additional SM processes, dark thermal masses and potential loops are **not included**. Thus we do not claim this EFT at arbitrary early-universe temperatures.

Check: an exact Z2 keeps phi stable, but permits phi-phi interactions and number-changing reactions; therefore the original closed-graph `N1+Nphi` can fail once new channels are added. The electron portal is elastic in chi1 number at tree level; other portals may violate that number. In chemical equilibrium a real self-conjugate phi has zero chemical potential when its number is not otherwise conserved.

## 1→2 decay: derived, not assigned

For chi1(P)→chi2(k)+phi(q), the parent-spin-averaged scalar Yukawa trace is

```text
|M|^2_avg = y^2 [ (m1 + m2)^2 - mphi^2 ].
lambda = [m1^2-(m2+mphi)^2] [m1^2-(m2-mphi)^2].
Gamma_rest = |M|^2_avg sqrt(lambda) / (16 pi m1^3).
```

The width vanishes at/below the physical threshold `m1 <= m2+mphi`, even at fixed nonzero y. For m2=mphi=0, it reduces to `Gamma=y^2*m1/(16*pi)`. The old BEST reference used a **prescribed** constant rest width: its convention is `sum |M|^2_BEST = g1 *16*pi*m1^3*Gamma/sqrt(lambda)`; using the trace above yields `g1 * y^2 [(m1+m2)^2-mphi^2]`, with `g1=4` for Dirac particles+antiparticles. The code tests this inverse mapping independently. **Beware:** this is a convention-matching exercise; channel multiplicity and degree-of-freedom mapping must still be checked in actual BEST source/rates before declaring new physical production.

Finite-T decay/inverse-decay is **not** obtained by multiplying the vacuum width by a guessed `gamma(T)`. The collision operator contains the complete gain–loss quantum factors

```text
f1 (1-f2) (1+fphi)  -  (1-f1) f2 fphi
```

and modified in-medium dispersion, phase space, and possible thermal masses. All terms must be inserted consistently. Momentum-dependent time-dilation alone is not a full thermal collision rate.

## Electron elastic scattering (tree level, massless e benchmark)

Take `s=(pchi+pe)^2`, `t=(pchi-pchi')^2`, `u=2*m1^2-s-t`, `D=s-m1^2`. Averaging the two spins of each incoming fermion and summing outgoing spins:

```text
|M|^2_avg = 2 G^2 [ D^2 + (u-m1^2)^2 + 2*m1^2*t ].
d sigma/dt = |M|^2_avg/(16*pi*D^2).
t_min = -D^2/s,   t_max=0.
sigma_total = G^2*D^2*(D^2+3*s^2)/(24*pi*s^3).
sigma_transport = int (1-cos theta_cm) d sigma
                = G^2*D^2*(3*D^2+4*s^2)/(48*pi*s^3).
```

The final two expressions are separately checked against numerical integration of the matrix element over t for several masses/energies. Near the non-relativistic chi rest limit `E_e << m1`, `sigma_transport -> G^2*E_e^2/(3*pi)`; at high energy `sigma_transport/sigma_total -> 7/8`. This is a perturbative tree-level cross section, with no electron mass, thermal screening, blocking, or mediator-propagator corrections.

## Rest-frame bath *transport collision frequency*, NOT kinetic relaxation

As a first analytic reference, consider one chi1 **at rest** in a massless, charge-symmetric e-/e+ bath, with `g_e=4`, zero chemical potential and initial-state FD occupations:

```text
nu_T(rest,T) = g_e/(2*pi^2) int_0^infinity dk k^2 fFD(k/T)
                                     sigma_transport(m1^2+2*m1*k)
nu_T(rest,T << m1) = 15*g_e*zeta(5)/(4*pi^3) * G^2 * T^5.
```

The integrator uses dimensionless `x=k/T`, and the asymptote has an independent analytic FD moment. This frequency **does not equal** `gamma_el(p,T)` entering the earlier number-conserving BGK operator: momentum/heat relaxation needs recoil, final-state Pauli blocking, a proper linearized Boltzmann/Fokker–Planck collision operator, distribution dependence, and physical multi-channel SM charge counting. Using `nu_T` as `gamma_el` now would be an unjustified model substitution. Electron mass, plasma effects and neutrino partners are absent. All exploratory ratios `nu_T/H` are labelled frequency ratios, **not** equilibration criteria.

## Independent control and accept/reject gates

1. Threshold test: width goes to zero with Källén phase space, unlike the previous fixed-`Gamma` mass-scan. Verify `Gamma=y^2*m1/(16*pi)` for massless daughters.
2. Convention test: independently invert BEST's existing `g1*16*pi*m1^3*Gamma/sqrt(lambda)` to the Yukawa width for massive daughters; do not rewrite historical BEST results.
3. Scattering test: compare analytical `sigma_total` and `sigma_transport` with independent quadrature over `d sigma/dt` for multiple `s,m,G`. Check positivity and angular limits.
4. Bath test: coupling-square scaling and heavy-rest FD asymptote; zero coupling must give zero.
5. Failure controls: unphysical thresholds, NaN and negative rates rejected; report all GeV dimensions.
6. **Next missing gates:** derive fully thermalized, Pauli-blocked `C_chi1[f]` from the same operator (including neutrinos where gauge invariance requires), compute number/energy moment and kinetic drag, compare to BGK with a calibrated `gamma`; add thermal masses and `g*_rho/g*_S`; independently validate BEST gamma/daughter distributions, long-time relic and observational claims. None of these are achieved by the current pilot.

Background motivation: [Yin, JHEP05(2023)180](https://arxiv.org/abs/2301.08735) and [Sakurai & Yin, JHEP03(2025)202](https://arxiv.org/abs/2410.18968). Existing BEST conventions: [best_massive_daughter_reference.py](best_massive_daughter_reference.py), [best_external_vacuum_reference.py](best_external_vacuum_reference.py). The vector-contact formula here is an explicit local tree-level trace/calculation, **not a published independent full plasma rate**.

## Reproducibility

```bash
python -m pytest -q tests/test_microphysics_portal.py
python kinetics/microphysics_portal.py
```

The command prints an **illustrative** benchmark (`m1=0.1 GeV, mphi=86e-9 GeV=86 eV`) with an arbitrary coupling, no assumption of observed ΩDM or freezeout. Its output is not an external BEST result; use independent CI or recomputation of this code before claiming numeric precision.


## Follow-up stage: thermal gain/loss and dimensional BEST bridge

[Microphysics thermal gain/loss](microphysics_thermal_gain_loss.py) checks the exact **on-shell occupation factors** from this Yukawa vertex, with explicit thermal Fermi and Bose functions. On an energy shell E1=E2+Ephi,

```text
D = f1(1-f2)(1+bphi), I=(1-f1)f2*bphi,
D / I = exp[(mu1-mu2-muphi)/T].
```

The test verifies `mu1=mu2` and `muphi=0` for equilibrium of a real self-conjugate scalar (if number changing is allowed), both signs of chemical affinity, strong low-energy Bose occupation and off-shell rejection. This is *not* a momentum-resolved thermal decay rate: matrix-element phase space, quantum statistical collision integrals and thermal masses remain necessary.

The same module exports `yukawa_best_bridge`: dimensionful GeV masses and the explicit Yukawa y become dimensionless `BEST masses=m/T`, `BEST Gamma=Gamma_rest/T`, and `BEST amplitude=|M|^2/T^2`, exactly inverted against the earlier independent analytic width. It deliberately does **not** execute BEST, relabel a fixed-width run as fixed-coupling, or insert the unblocked electron transport frequency into the BGK. A true pinned external BEST experiment is the next numerical gate after the microscopic collision channels and normalization are audited.
