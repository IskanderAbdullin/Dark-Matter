# Direct SPARC mass-model data audit: observational groundwork

**Date:** 2026-11-21. This is the first calculation on public observed SPARC velocity-table data in the continuation; no model fits or particle-mass posteriors are claimed.

## Source and exact schema
Official public source: https://astroweb.cwru.edu/SPARC/MassModels_Lelli2016c.mrt; Lelli, McGaugh & Schombert (2016), https://arxiv.org/abs/1606.09251. In this environment the 266-kB UTF-8 text was read through a public GitHub copy at https://github.com/LoaneShark/gravity-sandbox/blob/main/sparc/MassModels_Lelli2016c.mrt; Git blob SHA 3f0ea8e3577ec9e0226e4f6a65f1019d45f1aec6. Exact content agreement with the official server is not SHA256-verified. Cite original SPARC data producers, not the mirror.

The official text has per row: galaxy identifier, D [Mpc], R [kpc], Vobs [km/s], eVobs [km/s], Vgas, Vdisk, Vbul [km/s], and disk and bulge surface brightness. The reported eVobs omits certain inclination and distance systematic uncertainties.

## Computed results
A parsed text dataset of **3391 radii in 175 galaxies**, with 4 to 115 radii/galaxy (median 14), and **361 rows with signed negative Vgas**. Zero nonpositive radius, nonpositive velocity error or negative Vobs rows were detected.

For an *illustrative* stellar mass-to-light choice Udisk=0.5 and Ubulge=0.7, compute the nominal dark-halo velocity-squared residual:
\[
V_{\rm dm,nom}^2=V_{\rm obs}^2-V_{\rm gas}|V_{\rm gas}|
 -0.5 V_{\rm disk}^2-0.7V_{\rm bulge}^2. \tag{1}
\]
A negative Vgas means outward gas contribution at some radius: squaring it as a positive gas term without preserving the sign is a model error. There are **214** radii with negative nominal halo residual under these arbitrary fixed M/L values. This is not evidence of negative true dark-matter density: observational errors, uncertainty in M/L, inclination, distance and noncircular motions intervene. Do not delete these points or clip data values to zero; fit total predicted Vobs directly.

See sparc-qc-results.json for all 175 per-galaxy records and machine-checkable aggregate counts.

## What remains unverified
The 175 SPARC galaxies are not identical in measurements/scope to the original 22 LSB Halpha curves. The original sample's exclusion logic, all raw Halpha data, the SPARC quality flags and precise systematics need separate provenance work. At this stage there is NO fit to the Bose or Burkert models, NO common 86 eV posterior, NO reduced chi-square claim and NO cosmological constraint.

## Required statistical path
Obtain SPARC quality flags from Table1.mrt; model measured total velocities as
\[
V_{\rm model}=\sqrt{V_{\rm gas}|V_{\rm gas}|+
\Upsilon_dV_{\rm disk}^2+\Upsilon_bV_{\rm bul}^2+V_{\rm halo}^2},
\]
with physically valid positive total speed, galaxy-level correlated nuisance priors on M/L, distance and inclination, and actual measurement covariances where obtainable. For Bose, set nu from shared m,R,V0:
\[
\nu=\sqrt{\pi/2}\hbar^3/(G R^2 V_0 m^4), \quad 0<\nu\le\zeta(3/2).
\]
Compare with Burkert/PIS/NFW under the same likelihood and held-out predictive validation, not a ranking of mismatched reduced chi-square values. For lowered-Bose extension, use its own normalization and penalize extra cutoff degrees of freedom.
