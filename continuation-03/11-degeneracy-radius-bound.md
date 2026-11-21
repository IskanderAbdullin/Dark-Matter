# Algebraic degeneracy boundary as a galaxy-scale sanity check, not a mass measurement

In original JCAP 2021 for a single species with mass \(m\), central degeneracy satisfies
\[
\nu=n_0\lambda_{\rm th}^3=\frac{\sqrt{\pi/2}\hbar^3}{G\,R^2 V_0 m^4}
\leq\nu_c=\zeta(3/2).
\]
Use same convention \(V_0^2=kT/m\), \(R^2=V_0^2/(4\pi Gmn_0)\). For fixed \(m\) and \(V_0\), the purely **model-conditional** allowed scale radii satisfy
\[
R\ge R_{\min}(m,V_0)=\sqrt{\frac{\sqrt{\pi/2}\hbar^3}
{G V_0 m^4\nu_c}}\propto m^{-2}V_0^{-1/2}.
\]
For \(m=86\) eV and \(V_0=44.1\) km/s, with source units:
\[
\nu=\frac{41.8708571}{(R/\mathrm{kpc})^2(V_0/\mathrm{km/s})},\quad
R_{\min}\simeq 0.60286\;\mathrm{kpc}.
\]
Since the article rounds R=0.6 kpc and mass=86 eV, this proximity to the boundary is expected. Using the source's printed \(\nu=2.61\) in the exact relation gives \(m\simeq86.2246\) eV. This is an arithmetic consistency relation among independently fitted or assumed quantities, **not** a robust lower bound from observations.

Important: \(R\) is not necessarily a directly observed photometric/core radius; it is an internal equilibrium *scale parameter*. A physically truncated halo (lowered BE distribution) uses \(\nu_{\rm lowered}=4I(\psi_0;\eta)/\sqrt{\pi}\); the relevant cap differs and the shape relation between R and observed core radius changes. If \(V_0\) is mistaken for actual one-dimensional speed dispersion, the parameter inference may shift (true variance is \(V_0^2Li_{5/2}/Li_{3/2}\)). Do not apply the above inequality without rederiving it for a modified model.

Physical use: in a properly nuisance-marginalized population fit, report the fraction of galaxies whose best solutions lie on \(\nu_c\) and whether the apparent common mass is largely a **boundary-selection effect**. Compare profile likelihoods against priors for m, R, stellar mass/light and inclination, and test whether outliers require Bose-star parameters chosen after looking at residuals.
