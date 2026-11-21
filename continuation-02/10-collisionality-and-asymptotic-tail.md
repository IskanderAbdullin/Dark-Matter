# Note 10 — Collision rates and the outer asymptote of a critical noncondensed Bose halo
**Date 2026-11-21; numerical checks + conditional physical argument. Not a rigorous exclusion, not a real-galaxy refit.**

## 1. Why local thermodynamic equilibrium needs a mechanism
The original source assumes a Bose–Einstein distribution across the halo, including a chemical potential, while arguing contact self-interaction is negligible for the *equilibrium shape* of the noncondensed component. Negligible mean-field pressure and sufficiently rapid kinetic equilibration are distinct physical conditions. The observed virialization of a stellar or dark-matter system does not alone establish a fine-grained Bose–Einstein momentum distribution.

Take the source's published F583-4 \(R=0.6\) kpc, \(V_0=44.1\) km/s. For the source's cross-section scale \(\sigma/m=1.25\) cm²/g as an **illustrative upper-limit benchmark** (not a measurement), central density follows from its own scale equation:
\[
\rho_0=\frac{V_0^2}{4\pi G R^2}
=6.7648473\times10^{-24}\ {\rm g\,cm^{-3}}.
\]
The classical per-particle single-collision proxy is \(t_{\rm coll}(r)=1/\{\rho(r)(\sigma/m)V_0\}\), where \(V_0\) is used as an order-unity velocity proxy, rather than a properly integrated pairwise relative speed. The scattering mean free path \(\ell(r)=1/[\rho(r)(\sigma/m)]\) is \(38.32\) **kpc** at the center. Since \(R=0.6\) kpc, the central mean free path is ~64 times the scale radius. The corresponding \(t_{\rm coll}(0)=0.8497\) Gyr.

## 2. Density-dependent rates on the source-equation critical profile
The companion code solves the original ideal noncondensed spherical Bose–Poisson model at \(\nu=\zeta(3/2)\), normalized to the same \(R,V_0\); this is a close approximation to F583-4's tabulated rounded \(\nu=2.61\), but is NOT an independently fitted F583-4 density profile.

| \(r/R\) | \(\rho(r)/\rho_0\) | \(t_{\rm coll}\) [Gyr] | \(\ell\) [kpc] | Radial \(\tau_\sigma(0,r)\) |
|---:|---:|---:|---:|---:|
| 0.1 | 0.9463 | 0.898 | 40.50 | 0.00152 |
| 1 | 0.5855 | 1.451 | 65.46 | 0.01208 |
| 3 | 0.2286 | 3.717 | 167.63 | 0.02379 |
| 6 | 0.07602 | 11.178 | 504.14 | 0.03013 |
| 10 | 0.02655 | 32.007 | 1443.56 | 0.03298 |
| 20 | 0.005498 | 154.56 | 6971.10 | 0.03491 |

Here \(\tau_\sigma(0,r)=(\sigma/m)\int_0^r\rho(s)\,ds\) is an optical depth for **one radial passage**, NOT an accumulated orbital collision probability across a galactic lifetime. For any smaller \(\sigma/m\), mean free paths and times increase and optical depths decrease proportionally. Local occupation numbers can enhance Bose final-state scattering; collective gravitational encounters, condensation and transport require a proper kinetic treatment and could change these proxies. The collision time uses \(V_0\), not the Maxwellian mean two-body relative speed; therefore treat its numerical values as order-of-magnitude benchmarks, despite printed integration precision.

## 3. Separate asymptotic test of infinite untruncated mass
The ideal equilibrium equations are \((\xi^2 y')'=\xi^2 x\), \(x=\operatorname{Li}_{3/2}(z_0e^{-y})/\nu\). At large radii the Bose gas is dilute, so \(x\propto e^{-y}\). The asymptotic singular isothermal solution is
\[
x(\xi)\sim \frac2{\xi^2},\quad u^2(\xi)=\xi y'(\xi)\to2,
\quad M(<r)=\frac{V_0^2}{G}\,r\,u^2\sim\frac{2V_0^2}{G}r.
\]
The total mass diverges linearly, so a physical finite-radius halo must be matched to a tidal boundary or built from a truncated escape-energy distribution. Such truncation is a distinct *model*, not a numerical technicality.

I independently re-integrated the critical model from \(\xi=10^{-5}\) to \(\xi=10^5\), in \(t=\ln\xi\) using \(p=\xi y'=u^2\); then \(dy/dt=p\), \(dp/dt=e^{2t}x-p\). Outputs:

| \(\xi\) | \(\xi^2x\) | \(u^2\) |
|---:|---:|---:|
| 100 | 1.8257 | 1.9500 |
| 1,000 | 2.0570 | 2.0112 |
| 10,000 | 1.9817 | 1.9975 |
| 100,000 | 2.0051 | 2.0005 |

Small oscillatory approach is expected in the isothermal-envelope limit. The dimensions \(\xi=10^5\) are a mathematical numerical domain, not a physical halo extent.

## 4. Implications for papers 2 and 5
For paper 2: the untruncated equilibrium has infinite mass; explicitly implement a positive phase-space distribution with energy cutoff, a finite escaping population, and test stability.
For paper 5: show how the Bose gas acquires a local Bose–Einstein distribution on times much less than its age without arbitrarily selecting a large \(\sigma/m\) or assuming the condensate population; derive the nonlinear Boltzmann–Poisson or wave kinetic collision integral including Bose enhancement and compute orbital accumulation. A 0.03 single-pass optical depth alone cannot establish thermalization *or* exclude it.

## 5. Reproducibility and evidence
The companion Python file [collision_and_tail.py](code/collision_and_tail.py) depends on [mean_field_regularization.py](code/mean_field_regularization.py). It repeats density integration, collisional estimates, radial optical depths and log-radius tail; explicit numeric assertions check that the far field approaches \(2\) and that collision time at \(10R\) exceeds 13 Gyr for the fixed chosen proxy parameters. These tests passed locally. **No original observational velocities, hydrodynamic simulations, kinetic equation or self-consistent cosmic history were tested.**

Primary model: Abdullin & Popov 2021, https://doi.org/10.1088/1475-7516/2021/11/055. Relevant equilibrium caveat: Slepian & Goodman 2012, https://doi.org/10.1111/j.1365-2966.2012.21901.x. Later general review: Chavanis 2025, https://doi.org/10.3389/fspas.2025.1538434.