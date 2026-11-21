# Paper 4 — Population inference and falsification of noncondensed Bose halos with SPARC rotation curves
**Observational protocol manuscript, 2026-11-21. Not an executed fit.**
**Working title:** *Global particle-mass inference from baryon-aware rotation curves in noncondensed Bose–Einstein dark-matter halos*

## Abstract
We formulate a model-comparison framework capable of testing, rather than presupposing, the \(m\simeq86\) eV boson mass claimed from 20 selected low-surface-brightness galaxies by Abdullin & Popov (2021). The key improvements are a physically shared particle mass across the galaxy sample, explicit luminous-matter terms and nuisance uncertainties, a consistent rule for potential Bose-star populations, and galaxy-level out-of-sample tests. Modern SPARC observations, providing gas, disk and bulge contributions for 175 disk galaxies, are an appropriate public testbed. We derive a constrained reparameterization of the Bose halo model and specify posterior-predictive failures that would disfavor it even if a subset of galaxies has small reduced chi-squared.

## 1. Data provenance and source-model limitation
The source selected 22 LSB galaxies from the McGaugh et al. (2001) rotation-curve sample and neglected the baryonic component as relatively small. Of these, it identified two unsatisfactory galaxies (F571-8 and UGC-11648), omitted them from the common \(m\) scan that reported 86 eV, and improved their fits by adding a central population of condensed objects with extra model freedom. This makes the mass an important preliminary finding, not a bias-controlled population determination.

The [SPARC database](https://astroweb.cwru.edu/SPARC/) includes 175 late-type disk galaxies, high-quality H I/H-alpha rotation curves, gas contributions and 3.6 micron photometry with disk/bulge mass-model terms ([Lelli, McGaugh & Schombert 2016](https://arxiv.org/abs/1606.09251)). Eligibility cuts must be set without looking at fit residuals.

## 2. Physical halo parameterization
Let \(x_n(\xi;\nu)\) be the numerical solution of the noncondensed Poisson–Bose equation and \(u^2(\xi;\nu)=\xi^{-1}\int_0^\xi s^2x_n(s)\,ds\). Each galaxy has scale radius \(R_i\), velocity parameter \(V_{0,i}\), and degeneracy \(\nu_i\). The halo circular velocity is
\[
V_{{\rm h},i}(r)=V_{0,i}u(r/R_i;\nu_i),\quad
m=\left[\sqrt{\frac{\pi}{2}}\frac{\hbar^3}{G\nu_iR_i^2V_{0,i}}\right]^{1/4}. \tag{1}
\]
In the population model \(m\) is a *single global particle parameter*. Consequently the per-galaxy phase-space parameter is fixed by
\[
\nu_i=\sqrt{\frac\pi2}\frac{\hbar^3}{Gm^4R_i^2V_{0,i}},
\qquad 0<\nu_i\le\zeta(3/2). \tag{2}
\]
This is a nontrivial constraint coupling the galaxies; it must be enforced in the likelihood rather than fitting arbitrary \(\nu_i\) and reporting their implied masses afterward.

To add a Bose-star halo component, use
\[
x_{\rm BS}(\xi)=\alpha_i e^{-\gamma_i y_i(\xi)},\quad
y_i''+2y_i'/\xi=x_n+\alpha_i e^{-\gamma_i y_i}, \tag{3}
\]
or a physically more general distribution function. Unless a kinetic model predicts otherwise, \(\alpha_i\) and \(\gamma_i\) should **not** be set equal by fiat.

## 3. Correctly representing luminous matter
For radius \(r_{ij}\), SPARC decomposes the expected *squared* speed as
\[
V_{\rm mod}^2=V_{\rm halo}^2+\Upsilon_{d,i}V_{{\rm disk},ij}^2+
\Upsilon_{b,i}V_{{\rm bulge},ij}^2+
{\rm sgn}(V_{{\rm gas},ij})\,V_{{\rm gas},ij}^2. \tag{4}
\]
Some data files encode signed gas-velocity terms because the local radial gas contribution can be outward; squaring a signed gas velocity without its sign would produce a model error. The disk/bulge tables are normalized to a reference mass-to-light ratio and require proper scaling. Distances and inclinations change both the inferred radii and velocities; gas pressure support, asymmetric drift and non-circular motion matter for some low-mass systems.

A baseline point-level Gaussian likelihood is
\[
\log L=-\frac12({\bf V}_{\rm obs}-{\bf V}_{\rm mod})^\top C^{-1}
({\bf V}_{\rm obs}-{\bf V}_{\rm mod})-\frac12\log|2\pi C|. \tag{5}
\]
An actual analysis must propagate the catalog's measurement and astrophysical covariance, including possible inter-point correlations; (5) is a starting likelihood, not evidence that all uncertainties are independent.

## 4. Hierarchical inference without target selection
Construct \(\pi(m)\), astrophysically motivated priors on \(R_i,V_{0,i}\), mass-to-light ratios and kinematic nuisance variables, and a common *mixture or shrinkage* prior for whether each galaxy has a substantial condensed-object component. The binary/model selection should not label the source's two poor-fit galaxies as condensate-dominated before fitting. Separate terms for baryonic physics and profile-shape flexibility prevent an apparent “mass detection” from being a surrogate for unseen parameters.

A principled model hierarchy compares:
A. pure noncondensed Bose equilibrium;
B. the same plus one globally regularized Bose-star population model;
C. NFW with physically motivated baryonic response/feedback;
D. Burkert or pseudo-isothermal with matched priors;
E. an appropriate SIDM model, e.g. isothermal Jeans with baryons;
F. additional nonstandard models only if parameter accounting is fair.

Model preference is assessed by held-out galaxies, posterior predictive checks, stable evidence or predictive scores and residual morphology, not solely minimized \(\chi_{\rm r}^2\). Modern [2026 comparisons](https://doi.org/10.1093/mnras/stag2148) emphasize how distance and inclination uncertainties increase model degeneracy.

## 5. Identifiability and critical-boundary bias
The source table reports numerous best-fit \(\nu_i=2.61\), effectively the physical upper bound \(\zeta(3/2)\). Such boundary saturation requires profile-likelihood or posterior diagnostics, because a narrow best-fit mass can result from an imposed model boundary or informative priors. Likewise, as \(\nu\to0\), \(x_n\) approaches an isothermal sphere and mass ceases to be distinguishable from a classical shape. Show:
- complete marginalized and profile likelihoods for \(m\);
- fraction of galaxies with a boundary \(\nu_i\) posterior;
- sensitivity to LSB selection, baryonic priors and radial coverage;
- leave-one-galaxy-out predictive accuracy;
- injected synthetic galaxies with known \(m\) to test recovery and calibration;
- posterior predictive inner-slope and outer-radius diagnostics.

## 6. Falsifiable success and failure criteria
**Success requires** an astrophysically acceptable shared \(m\), residuals without systematic correlation with baryonic surface density, satisfactory predictions on held-out galaxies, and no unavoidable violation of independent cosmological constraints. Merely beating NFW on inner radii is insufficient because NFW is an incomplete baryon-affected comparison.

**Failure indicators** include a mass posterior highly sensitive to removing a few galaxies, systematic need for independently tuned \(\alpha_i,\gamma_i\) in most objects, observational preference for \(m\) incompatible with required cosmological streaming, or no improvement over simpler profiles at a fixed effective complexity.

## 7. Deliverables required before a scientific result
Archived sample definitions and raw data hashes; a public, tested solver; fitting scripts with likelihood/prediction checks; per-galaxy and aggregate posterior tables; blinded holdout comparisons; independent mock recoveries; code/data DOI; fully specified galaxy-quality cuts. **None of these fits has been run in this draft.** The paper specifies a falsifiable experiment and its likely degeneracies, not a newly measured boson mass.

Selected references: [Abdullin & Popov (2021)](https://doi.org/10.1088/1475-7516/2021/11/055); [SPARC data](https://astroweb.cwru.edu/SPARC/); [2026 SIDM Jeans alternative](https://doi.org/10.1093/mnras/stag969); [2026 galaxy-profile diversity](https://academic.oup.com/pasj/article-abstract/78/2/745/8514015).


## Appendix — first direct SPARC observed-data pilot (conditional and selection-biased)
A [new independent preliminary comparison](../continuation-02/13-sparc-20-galaxy-pilot.md) uses the publicly distributed SPARC mass-model table: 20 explicitly selected dark-matter-dominated galaxies, 353 fitting data points and 131 interleaved heldout radial points. The *original ideal noncondensed* Bose profile is constrained to a **single fixed mass m=86 eV** and two halo parameters (R,V0) per galaxy, enforcing the physical link \(\nu=\sqrt{\pi/2}\hbar^3/(GR^2V_0m^4)\le\zeta(3/2)\). Two-parameter Burkert and pseudo-isothermal (PIS) models use the same fiducial disk and bulge M/L (0.5 and 0.7) and the correct signed gas contribution. In this *specific* restricted fit, sum of heldout chi-square is Bose 616.83, Burkert 411.52 and PIS 252.95, so **the pseudo-isothermal model performs best**. This is NOT a proof of exclusion of an 86 eV boson: all 20 galaxies were selected through a nominal baryonic residual, M/L was fixed, distance/inclination systematic nuisance parameters and full covariances were not modeled, and the sample is not the 22 original JCAP Halpha points. The full SPARC sample with published quality metadata Q=1,2,3 and a priori inclination cuts must be fitted with honest hierarchical likelihoods before assessing an allowed mass range; see [quality metadata audit](../continuation-02/15-sparc-quality-and-sensitivity.md). A standalone Node reproduction script is included in the repository; its exact reproduction on another machine is an outstanding check.
