# Numerical note — free-streaming benchmark for a relativistically decoupled real Bose scalar

**Result: important quantitative alarm for the standard thermal version, NOT an observational bound.** Derived and evaluated 2026-11-21.

## Assumptions
One real scalar with occupation \(f_0(q)=1/(\exp(q/T_{b0})-1)\), zero primordial chemical potential, frozen comoving momentum after relativistic decoupling. A single boson mass \(m\), no post-decoupling condensate, no elastic collision term in subsequent propagation, no decay, no entropy injection. For 86 eV use \(g_{*S,d}=55\), \(T_{\nu0}=1.95\) K, \(T_{b0}=T_{\nu0}(10.75/g_{*S,d})^{1/3}=9.7519\times10^{-5}\) eV (energy units).

**Background:** flat FLRW with \(H_0=67.4\) km/s/Mpc, \(\Omega_m=0.315\), \(\Omega_r=9.2\times10^{-5}\), \(\Omega_\Lambda=1-\Omega_m-\Omega_r\); held fixed for the calculation. A fully self-consistent background should adjust radiation for the modeled scalar.

## Formula evaluated — not just speed at equality
With dimensionless comoving momentum \(x=q/T_{b0}\), the number-weighted mean physical speed (units of \(c\)) is
\[
\langle\beta(a)\rangle_n =
\frac{\int_0^\infty dx\,\frac{x^2}{e^x-1}\,
    \frac{xT_{b0}}{\sqrt{(xT_{b0})^2+(am)^2}}}
{\int_0^\infty dx\,\frac{x^2}{e^x-1}}.
\]
The comoving ballistic distance traveled by \(a\) is
\[
\lambda_{\rm ballistic}(a)=\int_{10^{-12}}^{a}
\frac{c\langle\beta(a')\rangle_n}
{a'^2 H(a')}\,da'
=\frac{c}{H_0}\int_{10^{-12}}^a
\frac{\langle\beta(a')\rangle_n\, da'}
{\sqrt{\Omega_r+\Omega_m a'+\Omega_\Lambda a'^4}} .
\]
Its use of a mean velocity produces a *representative ballistic distance*, not an exact phase-space-response kernel; no transfer function is inferred from it. The momentum integral is normalized by \(2\zeta(3)\).

## Calculated numbers for 86 eV
- \(\langle q\rangle/T_{b0}=2.70117803\); characteristic nonrelativistic transition \(z_{\rm nr}\approx326,478\).
- \(\Omega_bh^2=0.120314\); \(\Delta N_{\rm eff}=0.064817\) in the one-real-scalar entropy-conservation benchmark.
- \(z=3400\) (equality): \(\langle v\rangle_n=3122.57\) km/s; \(\lambda_{\rm ballistic}=6.6583\) comoving Mpc.
- \(z=1100\) (recombination): \(1011.00\) km/s; \(7.6202\) comoving Mpc.
- \(z=8\): \(8.2645\) km/s; \(9.0097\) comoving Mpc.
- \(z=5\): \(5.5097\) km/s; \(9.0364\) comoving Mpc.
- \(z=0\): \(0.9183\) km/s; \(9.1001\) comoving Mpc.

The code was tested with a coarser momentum/scale-factor grid. Difference in the distance to equality was below 0.5%; analytic BE mean momentum agreed within \(5\times10^{-5}\). This is a numerical integration cross-check, not independent data validation.

## Same assumed present density, several masses
For each row, choose \(g_{*S,d}=55(m/86\,\mathrm{eV})\) so the **simplified** thermal abundance is held at \(\Omega_bh^2\simeq0.120314\). The implied decoupling temperatures need a separate particle-physics model and are not all automatically viable.

| m [eV] | g*S,d | ballistic distance by equality [comoving Mpc] | by z=5 [comoving Mpc] |
|---:|---:|---:|---:|
| 29 | 18.55 | 19.73 | 29.86 |
| 50 | 31.98 | 11.63 | 16.53 |
| 86 | 55.00 | 6.66 | 9.04 |
| 152 | 97.21 | 3.62 | 4.73 |

## What this can and cannot establish
It demonstrates that an **ordinary thermal noncondensed relic** in the mass range implicated by the JCAP source has long ballistic propagation paths in a standard cosmological background. It does **not** derive the suppression of the linear matter power spectrum (which depends on the full Boltzmann hierarchy, gravitational potentials, collisions, and the distribution's shape). It cannot directly be compared numerically with a thermal-*fermion* WDM 5.7 keV limit without calculating the corresponding Bose transfer function and likelihood. Moreover, it is unrelated to early cold-condensate or nonthermal scenarios unless their actual \(f(q,a)\) is specified.

**Next reproducible requirement:** implement the exact model in CLASS/CAMB (or independent collisionless Einstein–Boltzmann hierarchy for frozen \(f_0(q)\)), compute \(T(k)\), \(k_{1/2}\), CMB and Lyman-alpha flux-power predictions and vary the primordial cold fraction/time of formation; evaluate \(\Omega h^2\) and \(\Delta N_{\rm eff}\) simultaneously. An assumed galactic BE equilibrium at late times is not a substitute for this primordial calculation.

Sources: Abdullin & Popov 2021 https://doi.org/10.1088/1475-7516/2021/11/055; Rogers et al. (2024) https://doi.org/10.1103/PhysRevD.109.043511; mixed cold+warm Lyman-alpha 2025 https://doi.org/10.1103/4k29-h99l; SPARC https://arxiv.org/abs/1606.09251.