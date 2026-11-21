# Paper 3 — Cosmological consistency of 50–150 eV bosons with galaxy-scale Bose-degenerate halos
**Working theory paper, 2026-11-21; no calculated cosmological likelihood or allowed parameter window.**
**Working title:** *Can eV–hectoelectronvolt bosonic relics be cosmologically viable when a noncondensed halo component dominates?*

## Abstract
An 86 eV mass selected from dark-matter-only low-surface-brightness galaxy fits does not, by itself, demonstrate that such bosons could form the observed cosmic density perturbations. We derive a normalized thermal-relativistic-relic abundance benchmark for a single real scalar, distinguish that benchmark from cold-condensate and nonthermal production, and define a common phase-space-to-linear-perturbations test for the original two-component Bose-halo proposal. The controlling prediction is not merely the inner rotation-curve shape but the joint relic density, radiation contribution, particle velocity distribution, linear transfer function, and cosmological timing of Bose-star formation. A strong small-scale-structure tension exists for a freely streaming 86 eV thermal population; its severity for Bose distributions or mixed early condensed phases must be computed rather than borrowed directly from fermionic thermal warm-dark-matter limits.

## 1. Why the source mass is not a primordial mass determination
Abdullin and Popov (JCAP 2021 055) use a Thomas–Fermi Bose–Poisson halo with a noncondensed central degeneracy parameter \(\nu\le\zeta(3/2)\), plus small Bose condensates grouped into compact objects. Its joint mass scan over a selected 20 of 22 LSB galaxies reports \(m=86\) eV. This is an inference in a stationary galaxy model: it neither fixes the primordial distribution \(f(q,a)\) nor proves that the required chemical potential and macroscopic temperature arise dynamically. In particular, “cold” condensation at late times cannot undo the suppression of primordial density perturbations that has already occurred through early free streaming.

## 2. A controlled abundance benchmark
For **one real scalar degree of freedom**, relativistically decoupled with negligible chemical potential, with no macroscopic condensate and with conserved Standard Model comoving entropy,
\[
T_b/T_\nu=(10.75/g_{*S,d})^{1/3}. \tag{1}
\]
Because \(n_b=\zeta(3)T_b^3/\pi^2\) while one neutrino flavor plus antineutrino has \(n_{\nu+\bar\nu}=(3/2)\zeta(3)T_\nu^3/\pi^2\), the standard neutrino reference formula implies
\[
\Omega_b h^2\simeq \frac23\frac{m}{93.14\ \mathrm{eV}}\frac{10.75}{g_{*S,d}},\qquad
\Delta N_{\rm eff}=\frac47(T_b/T_\nu)^4. \tag{2}
\]
These benchmark expressions **do not apply unchanged** to a complex scalar with conserved asymmetry, post-decoupling entropy injection, a nonzero chemical potential, cannibal processes or an early BEC.

For the source's representative mass \(m=86\) eV, (2) gives \(g_{*S,d}=25 \Rightarrow \Omega h^2\simeq0.265\) and \(\Delta N_{\rm eff}\simeq0.185\); \(g_{*S,d}\simeq55 \Rightarrow\Omega h^2\simeq0.120\) and \(\Delta N_{\rm eff}\simeq0.065\). The exact decoupling epoch and entropy degrees of freedom must be consistent with the particle interactions and number of components; these values are arithmetic checks, **not** a fit to CMB/BBN data.

## 3. Perturbation evolution as a decisive test
For a frozen collisionless comoving momentum \(q=ap\), the physical speed of a particle is
\[
v(q,a)=\frac{q}{\sqrt{q^2+a^2m^2}}
\]
in natural units. For a distribution \(f_0(q)\), characterize its physical free-streaming history by the appropriately momentum-weighted integral
\[
\lambda_{\rm fs}(a)\sim\int_{a_d}^a\frac{\langle v(q,a')\rangle_{f_0}}{a'^2H(a')}\,da'. \tag{3}
\]
The full linear transfer function, including metric perturbations, time-varying sound speeds, and any Bose-fluid/condensate component, requires a Boltzmann solver such as CLASS or CAMB modified for the **actual** \(f_0(q)\). Formula (3) is a diagnostic, not a substitute for cosmological perturbation calculations.

If a fraction \(f_c(a)\) becomes cold early, its gravitational coupling with the warm component changes the total growth and the total \(T(k)\). If the same fraction condenses only after the warm component erases small-scale modes, it cannot regenerate a generic CDM-like initial halo mass function without an additional nonlinear mechanism. Such timing is the central survival requirement for the original model.

## 4. Modern external evidence and scope of constraints
A 2024 Lyman-alpha analysis reports an approximate 95% confidence lower limit of **5.7 keV** on a specified *thermal relic WDM* template; under narrower matched priors it reports **4.1 keV**. This is [Phys. Rev. D 109, 043511 (2024)](https://doi.org/10.1103/PhysRevD.109.043511). Mixed warm-plus-cold models have distinct allowed fractions, as demonstrated by [2025 high-redshift Lyman-alpha constraints](https://arxiv.org/abs/2504.06367). For distribution-specific modeling see [nonthermal WDM limits (2023)](https://doi.org/10.1103/PhysRevD.108.043518). It would be invalid to assert simply “\(m_b>5.7\) keV” by copying a fermionic thermal template bound onto a Bose gas with a nonthermal or condensed component. It would be equally invalid to conclude that quantum statistics alone exempts an 86 eV relic from free streaming.

## 5. Four mutually distinct model scenarios
**T0: Standard relativistically decoupled real scalar, no condensate.** Determine its abundance, \(N_{\rm eff}\), linear matter transfer function and Lyman-alpha flux. This is the cleanest null benchmark.

**T1: Primordial asymmetry / early cold condensate.** Solve the finite-temperature chemical potential and charge-density evolution; identify the condensate formation epoch relative to horizon entry of dwarf-scale modes. Show an explicit particle interaction or conserved-number mechanism; otherwise this is not a complete theory.

**T2: Nonthermal cold phase-space production.** Specify decay/freeze-in production and construct its *normalized* \(f(q)\), then derive its growth predictions. The mass 86 eV alone does not set streaming: momentum matters.

**T3: Mixed cosmology and late Bose-star nucleation.** Compute kinetic relaxation and early dark-matter perturbations together. Late compact-object formation can alter halo cores and satellite observables but does not automatically repair an erased primordial power spectrum.

These scenarios must not be interchanged when interpreting an observed \(m\) posterior.

## 6. Implementable inference and stop conditions
1. Fix the microscopic degrees of freedom, chemical potential, interactions and production channel, and require \(\Omega_{\rm DM}h^2\) to match cosmological data within model uncertainties.
2. Compute \(f_0(q)\), \(\Delta N_{\rm eff}\), and CMB/BBN consistency (including a possible dark radiation bath).
3. Integrate linear Boltzmann perturbations to obtain \(T(k)\) and \(P(k,z)\) at Lyman-alpha redshifts; propagate IGM uncertainties rather than reading a single published mass limit as universal.
4. Compare halo mass functions, dwarf counts, reionization and lensing, while tracking a time-dependent condensate fraction.
5. Only then insert the surviving mass/coupling combinations into baryon-aware galactic equilibrium fits.
6. **Stop condition:** if no dynamically derived \(f_0(q)\) and condensate history gives acceptable small-scale \(P(k)\) at \(m\sim50\)–150 eV, the original thermal interpretation must be rejected or altered, even if rotation curves fit.

## 7. What has and has not been established
Equations (1)–(3) and the numerical abundance examples are derived benchmark results. A credible, observationally allowed 86 eV boson model has **not** been demonstrated. The cosmological likelihood, transfer function, detailed thermal history, particle couplings and condensation epoch remain open calculations. Relevant context: [Boyanovsky et al. (2008)](https://doi.org/10.1103/PhysRevD.77.043518); [Abdullin & Popov (2021)](https://doi.org/10.1088/1475-7516/2021/11/055); [Chavanis (2025) review](https://doi.org/10.3389/fspas.2025.1538434).


## Appendix — a thermal-relic velocity benchmark at matter–radiation equality (new arithmetic check)
A simple numerical scale makes the cosmological viability question concrete. Take exactly the *single relativistically decoupled real scalar benchmark* of Section 2, with zero primordial chemical potential, \(m=86\,{\rm eV}\), \(g_{*S,d}=55\), present neutrino temperature \(T_{\nu0}=1.95\,{\rm K}\), and approximate equality redshift \(z_{\rm eq}=3400\). For a massless Bose–Einstein momentum distribution prior to kinetic decoupling the **mean physical momentum** is \(\langle p\rangle=\pi^4/[30\zeta(3)]\,T_b\simeq2.70118\,T_b\), with \(k_B=1\) in this paragraph. From entropy conservation,
\[
T_{b0}\simeq9.7519\times10^{-5}\ {\rm eV},
\quad
\langle p\rangle_{\rm eq}\simeq 2.70118\,T_{b0}(1+z_{\rm eq})
\simeq0.896\ {\rm eV}.
\]
Therefore \(\langle p\rangle_{\rm eq}/\sqrt{m^2+\langle p\rangle_{\rm eq}^2}\simeq0.01042\), or roughly **\(3.12\times10^3\ {\rm km\,s^{-1}}\)**. This speed is computed from **mean momentum**; it is only a rough proxy for the mean speed, but is accurate at leading nonrelativistic order \(p/m\ll1\). It is not a 1D velocity dispersion and must not be substituted directly into a Boltzmann transfer function. \(m/(\langle p\rangle_0)\simeq 3.26\times10^5\) is the approximate redshift at which the typical momentum first becomes nonrelativistic. The same benchmark's \(\Omega h^2\simeq0.1203\) and \(\Delta N_{\rm eff}\simeq0.0648\) are internally consistent with Section 2 arithmetic.

**Inference:** even a relic whose abundance happens to match the observed dark-matter abundance at this decoupling temperature can retain appreciable thermal velocities near equality. This strengthens the case for explicit small-scale transfer-function calculations. The result applies **only to the stated thermal relic**; early condensates or nonthermal low-momentum spectra can behave differently and must be derived microscopically.


## Appendix B — full ballistic streaming integral for one specified 86 eV thermal real scalar
The [continuation numerical calculation](../continuation-02/09-thermal-streaming-benchmark.md) carries out the number-weighted momentum integral for a frozen relativistic Bose–Einstein distribution and integrates the resulting mean particle speed through an explicitly specified flat expanding cosmology. **This is not a Boltzmann transfer function.** In the scenario \(m=86\) eV, \(g_{*S,\rm dec}=55\), \(H_0=67.4\) km/s/Mpc, \(\Omega_m=0.315\), \(\Omega_r=9.2\times10^{-5}\), it yields number-averaged \(v(z=3400)\simeq3123\) km/s and a **ballistic** comoving traveled distance about **6.66 Mpc** up to equality, **9.04 Mpc** up to \(z=5\), **9.10 Mpc** to today. At this decoupling ratio \(\Omega_b h^2\simeq0.1203\), \(\Delta N_{\rm eff}\simeq0.0648\), with the background expansion held fixed rather than refitted to extra radiation. Distances are not direct half-mode scales or Lyman-alpha exclusions. An early cold condensate, a nonthermal spectrum or substantial scattering require a completely different calculation. These results should be treated as a first-priority *consistency alarm*, not a quantitative observational rejection.
