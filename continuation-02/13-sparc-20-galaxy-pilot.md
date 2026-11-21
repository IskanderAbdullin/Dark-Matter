# Preliminary 20-galaxy SPARC pilot: 86 eV Bose vs Burkert vs PIS

Date 2026-11-21. **A first actual data fit in an explicitly simplified model**, not a publication-grade constraint or universal rejection of the original theory.

## Source and selection
Direct public SPARC mass model measurements (Lelli, McGaugh & Schombert 2016), https://astroweb.cwru.edu/SPARC/MassModels_Lelli2016c.mrt, loaded in this environment from independently accessible public mirror https://github.com/LoaneShark/gravity-sandbox/blob/main/sparc/MassModels_Lelli2016c.mrt, source blob SHA 3f0ea8e3577ec9e0226e4f6a65f1019d45f1aec6. 3,391 points, 175 galaxies overall. For **this pilot only**, deterministically sort galaxies alphabetically and choose first 20 with (i) at least 15 radii, (ii) median nominal dark velocity-squared fraction >0.7, and (iii) nonnegative nominal residual at every radius. Selection is **data-dependent**, so results are not an unbiased test of the population.

20 selected galaxies, 484 measured radii: 353 fit points, 131 holdout points (every fourth sorted radius, beginning at index 0). Split is radial interleaving, not an independent galaxy sample; neighboring observations may have correlated residuals.

## Three physically defined shapes, two fitted parameters each
1. Original *untruncated ideal* noncondensed Bose profile, constrained to single globally **fixed, not fitted**, particle mass m=86 eV. Free per-galaxy R in kpc and V0 in km/s. Correct degeneracy
\[
\nu=\frac{\sqrt{\pi/2}\hbar^3}{G R^2 V_0 m^4}
   =\frac{41.8708571114}{[R/{\rm kpc}]^2[V_0/{\rm km\,s^{-1}}]},\ 0<\nu\le\zeta(3/2).
\]
For each nu the dimensionless Poisson–Bose ODE is solved in log radius, 102 nu templates and interpolated u²(xi)=xi*y'. The near-singular polylog is computed by series/asymptotic; independent fixed-profile checks versus previously calculated values at nu=0.03 and critical nu passed.
2. Burkert density: two free scale parameters. Circular speed shape proportional to [2ln(1+x)+ln(1+x²)−2atan(x)]/x.
3. Pseudo-isothermal: two free parameters, speed shape 1−atan(x)/x.

For ALL models disk M/L =0.5, bulge M/L=0.7 fixed throughout, gas sign kept: V_pred²=Vgas|Vgas|+0.5 Vdisk²+0.7 Vbul²+V_halo². Parameter sets with V_pred²<0 at ANY observed radial position were explicitly prohibited (including heldout radii, but without using the heldout observed speeds). Total predicted V and random eVobs were compared by the diagonal weighted chi-square. All optimizations used coarse log grids then three-start Nelder–Mead. No M/L priors, covariances, distance/inclination, pressure-support treatment, extra Bose-star population or special selection penalties.

## Output and interpretation
| Model | Train chi² / 353 pts | Holdout chi² / 131 pts | Best holdout fit in galaxies |
|---|---:|---:|---:|
| Bose, fixed 86 eV | 1749.95 | 616.83 | 5/20 |
| Burkert | 961.84 | 411.52 | 4/20 |
| Pseudo-isothermal | **550.79** | **252.95** | **11/20** |

These are sums of **ordinary chi-square residuals**, not sums of published reduced chi-squares. They do not incorporate unknown observational covariances or additional systematic errors; many results fit poorly in absolute terms. The **conditional** result is that PIS describes/predicts these selected observed curves better under the deliberately restricted fitting assumptions. It is NOT a likelihood-derived mass exclusion, nor proof that a finite-radius Bose profile or a two-component Bose-star halo fails.

The previous unphysical holdout issue was fixed before reporting these values: models must be real/positive on every radial point. Source analysis in code/sparc_pilot_compare.js reproduces the pipeline; it reads the official text file from a local path and outputs all per-galaxy parameters and holdout scores.

## Major required robustness checks
- Download official data from origin and compare cryptographic hashes against public mirror; add original per-galaxy quality and source metadata.
- Fit with stellar M/L nuisance parameters including priors (and concentration/stellar population uncertainties), and distance/inclination with their error models.
- Extend to all 175 galaxies, not this selected 20, and separately reproduce the ORIGINAL 22-galaxy JCAP sample.
- Use full covariance/extra scatter, posterior predictive checks, stratified holdout by galaxy, and profile likelihoods for global m (rather than fixing m=86).
- Repeat with phenomenological Bose-star component, plus the newly solved lowered-Bose finite-mass solution; pay for added parameters/priors.
- Interpret the cosmological Boltzmann transfer and kinetic relaxation conditions before claiming dark-matter viability.
