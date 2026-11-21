# Paper 5 — Kinetic production and gravitational tests of compact Bose-star populations
**Theory/observational strategy draft, 2026-11-21. No validated mass function or detection forecast.**
**Working title:** *Population kinetics, equilibrium backreaction and extended-lens tests of Bose stars in 50–150 eV bosonic dark-matter halos*

## Abstract
The condensed component of a bosonic dark-matter halo affects rotation curves only through its mass-weighted spatial distribution, yet the underlying individual Bose-star masses, radii, abundances, production rates and survival probabilities are unknown in the 2021 model. We promote that unobserved population from an empirical density term to a kinetic mass function, establish conditions for contact-relaxation relevance, and derive joint lensing, dynamical-heating and host-halo observables. We show explicitly why a radius bound alone cannot predict microlensing, and formulate falsifiable links among condensate formation, spatial segregation and observational limits.

## 1. Source hypothesis and missing scales
Abdullin & Popov (2021) propose that noncondensed \(\mathcal{O}(10^2\,{\rm eV})\) bosons dominate much of a galaxy halo, while compact condensates (“Bose stars”) can cluster toward the center and correct the inner rotation curves of otherwise poorly fit galaxies. Their repulsive self-interaction Thomas–Fermi radius bound is
\[
R_{\rm BS}\lesssim83.7\,{\rm AU}
\left[\frac{(\sigma/m)}{1.25\ {\rm cm^2\,g^{-1}}}\right]^{1/4}
\left(\frac{m}{50\ {\rm eV}}\right)^{-5/4}. \tag{1}
\]
This bound is specific to a particle/cross-section/approximation regime. It is not an object mass prediction, abundance estimate, formation time, stability criterion, or lensing exclusion. A non-interacting star follows a different mass–radius relation and may be much smaller.

## 2. Population state variable
Let \(\mathcal{F}(M,r,t)=d^2N_{\rm BS}/(dV\,dM)\) be a number density per mass interval (possibly velocity-resolved). Then
\[
\rho_{\rm BS}(r,t)=\int M\mathcal{F}(M,r,t)\,dM,\quad
f_{\rm BS}(r,t)=\frac{\rho_{\rm BS}}{\rho_{\rm BS}+\rho_n}. \tag{2}
\]
The central phenomenological parameters of the source are \(\alpha=\rho_{\rm BS}(0)/\rho_n(0)\) and \(\gamma=V_0^2/v_{\rm BS}^2\). They are separate dynamical quantities. The source's \(\alpha\approx\gamma\) is a fitting simplification, not a consequence of the Jeans equation alone. A successful population model should *predict* their relation and galaxy-to-galaxy scatter.

An effective mass-function kinetic equation is schematically
\[
\frac{\partial \mathcal{F}}{\partial t}
+\nabla_r\!\cdot(\mathbf v\mathcal{F})
-\nabla_r\Phi\cdot\nabla_v\mathcal{F}
=C_{\rm nucleation}+C_{\rm coagulation}-C_{\rm disruption}-C_{\rm evaporation}. \tag{3}
\]
Each collision operator needs a specified particle interaction, condensate structure and dark-matter phase-space distribution. Conservation of boson number, mass and energy between the smooth gas and bound stars must be checked explicitly.

## 3. A physical relaxation clock
For contact scattering with density \(\rho\), speed \(v\) and scattering cross-section per particle mass \(\sigma/m\),
\[
\Gamma_{\rm coll}\sim \rho(\sigma/m)v,\qquad
t_{\rm coll}\sim [\rho(\sigma/m)v]^{-1}. \tag{4}
\]
At \(\rho=10^{-24}\ {\rm g\,cm^{-3}}\), \(v=100\ {\rm km\,s^{-1}}\) and \(\sigma/m=1.25\ {\rm cm^2\,g^{-1}}\), this crude binary collision time is **2.54 Gyr**. It becomes longer for a smaller cross-section or lower density. This is neither a gravitational two-body relaxation time nor the Bose condensation time: occupancy factors, binding, many-body scattering and wave effects may modify the kinetic equation, and no numerical condensation epoch follows from (4). The 2023–2024 kinetic studies [multicomponent Bose stars](https://doi.org/10.1103/PhysRevD.108.043535) and [self-interacting nucleation](https://doi.org/10.1103/PhysRevD.109.016002) motivate careful treatment of applicable mass and occupation regimes; their fuzzy-dark-matter results cannot simply be copied to 86 eV.

## 4. Gravitational lensing and finite-object size
For a lens object with mass \(M\) at distance \(D_L\) and background source at \(D_S\), the characteristic point-mass Einstein radius in the lens plane is
\[
R_E=\sqrt{\frac{4GM}{c^2}\frac{D_LD_{LS}}{D_S}}. \tag{5}
\]
A boson star behaves approximately as a point microlens only if its projected mass distribution is compact compared with the relevant \(R_E\) and source size. If \(R_{\rm BS}\gtrsim R_E\), extended-lens light curves and astrometric deflections are needed. Therefore (1) does not tell us whether microlensing surveys see the population: the crucial missing variable is \(\mathcal F(M,r,t)\).

A model event-rate calculation must integrate the physical extended-lens cross-section over \(\mathcal{F}(M,r,t)\), transverse velocity and lens–source geometry, then fold in survey cadence and selection efficiency. Consult [extended-object microlensing (2024)](https://doi.org/10.1103/PhysRevD.109.123004), [Gaia astrometric extended-lens study (2026)](https://doi.org/10.1088/1475-7516/2026/09/072), and [axion-star mass-function study (2024)](https://doi.org/10.1088/1475-7516/2024/08/023). These sources establish relevant methods, not direct mass-86-eV exclusions.

## 5. Independent dynamical probes
A compact-object subpopulation can produce stochastic gravitational forces and tidal perturbations affecting wide stellar binaries, streams or disks. Observables depend on the combination of object mass and number density, not solely the smooth summed density \(\rho_{\rm BS}\). Use the [wide-binary minicluster constraints (2025)](https://doi.org/10.1103/PhysRevD.111.043042) as a methodology reference, with object-specific velocity and extended-mass effects rederived. Also test central mass clustering against baryon-aware rotation curves, stellar heating and strong-lens flux anomalies; avoid double counting correlated observations.

The proposed population can generate a steep effective inner halo profile through mass segregation, but a stationary isotropic Jeans fit does not show that such segregation *forms* within a Hubble time. A time-dependent multi-component gravitating simulation is therefore essential.

## 6. Numerical programme
1. Map allowed boson-star equilibrium branches across \(m,\sigma/m,M\) with Schrödinger–Poisson or Gross–Pitaevskii–Poisson as appropriate; impose relativistic stability when needed.
2. Calculate nucleation/relaxation rates in the **noncondensed 50–150 eV** regime using relevant occupation numbers and cross-sections; do not extrapolate fuzzy-ultralight scaling blindly.
3. Evolve the mass function while conserving mass and energy; calibrate evaporation and merger prescriptions on resolved simulations.
4. Evolve \(\rho_n+\rho_{\rm BS}+\rho_b\) self-consistently; predict \(\alpha,\gamma\) and radial condensate fraction, *not* fit them independently per galaxy by construction.
5. Simulate finite-size lensing event rates and stellar dynamical effects; compare with likelihoods from Gaia/OGLE/HSC and other surveys where the sensitivity is actually nonzero.
6. Feed all viable population models back into the cosmological transfer-function audit and SPARC galaxy likelihood.

## 7. Discriminating outcomes and uncertainty ledger
**Positive prediction (conditional):** a real central condensate population large enough to affect galaxy rotation curves also supplies a nonzero mass-function-weighted signal in lensing or dynamical perturbations; whether accessible depends on \(\mathcal F\), object compactness and exposure.

**Falsification paths:** a population required by kinematics violates quantitative extended-lens or wide-binary constraints; condensates cannot form or survive by the measured galaxy age; required \(\alpha,\gamma\) vary without plausible dynamics; microphysics needed for condensation violates self-interaction bounds; or cosmology excludes the underlying primordial distribution.

**Not yet computed:** Bose-star mass function, formation epoch, observable event rate, allowed cross-section range, and decisive observational upper limits. Those cannot responsibly be replaced by invented numerical findings. Prior work: [Abdullin & Popov (2021)](https://doi.org/10.1088/1475-7516/2021/11/055); [2025 BECDM review](https://doi.org/10.3389/fspas.2025.1538434); [2024 axion-star plus PBH lensing](https://doi.org/10.1088/1475-7516/2024/10/013).


## Appendix — collision-rate and optical-depth diagnostic of Bose-halo equilibrium
A [separate profile calculation](../continuation-02/10-collisionality-and-asymptotic-tail.md) finds, for the source galaxy F583-4's \(R=0.6\) kpc, \(V_0=44.1\) km/s, and a deliberately large illustrative \(\sigma/m=1.25\,{\rm cm}^2/{\rm g}\): central \(\rho_0=6.765\times10^{-24}\,{\rm g\,cm^{-3}}\), classical central mean free path \(\ell_0=38.3\) kpc and collision-time proxy \(t_0=\ell_0/V_0\simeq0.85\) Gyr. On the ideal critical profile, \(t_{\rm coll}\sim11.2\) Gyr at \(6R\), \(32.0\) Gyr at \(10R\); the radial single-pass optical depth from the center to \(10R\) is just \(0.033\). These are NOT complete Bose-enhanced quantum kinetic rates and do not by themselves rule out a common equilibrium distribution, but make the dependence on microscopic thermalization and orbital histories impossible to omit. The use of the 1-D \(V_0\) in place of an average pair relative speed adds another order-unity rate uncertainty. Contact cross section here is a chosen upper-limit-scale proxy, not a measurement or robust prediction of the original model.
