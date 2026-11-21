# Paper 2 extension — a phase-space-positive lowered Bose–Einstein halo of finite radius and mass
2026-11-21; independently derived ansatz and numerically integrated equilibrium. **Not** a proof of kinetic realization, stability, galaxy fits or novel statistical mechanics.

## 1. Motivation
The original JCAP model uses a Bose–Einstein distribution depending on the full energy E=p²/(2m)+mPhi. Its dilute outskirts approach the isothermal r^-2 tail (verified in note 10), so total mass grows without bound. Galaxies should instead be embedded in a finite escape/tidal environment. A spatial density taper need not correspond to stationary, nonnegative phase space.

## 2. A positive energy-cutoff distribution
For cutoff E_t and global chemical-potential parameter mu, define
\[
 f(E)=\begin{cases}
 (e^{(E-\mu)/k_BT}-1)^{-1}-(e^{(E_t-\mu)/k_BT}-1)^{-1},&E<E_t\\
 0,&E\ge E_t .
 \end{cases} \tag{1}
\]
Up to the usual quantum phase-space cell normalization, this is a King-like **lowered Bose** function. It is nonnegative when mu<E<E_t; the phase-space population vanishes continuously at escape. Because it depends only on conserved single-particle E, it solves the collisionless stationary Vlasov equation in a time-independent potential. This does NOT imply strict collisional Bose–Einstein equilibrium: subtraction of escape occupation is a phenomenological ansatz.

Define V0²=kBT/m, psi=(E_t-mPhi)/kBT, eta=(E_t-mu)/kBT, q=p/sqrt(2mkBT). Require 0<psi0<eta, with
\[
 f(q;\psi,\eta)=\frac1{e^{\eta-\psi+q^2}-1}-\frac1{e^\eta-1},
 \quad 0\le q\le\sqrt\psi;\quad f=0 \text{ otherwise}. \tag{2}
\]
Then 0<=q²<=psi ensures the first denominator <= the second, so f>=0.
Density is proportional to
\[
 I(\psi;\eta)=\int_0^{\sqrt\psi}q^2f(q;\psi,\eta)\,dq,\qquad x(\psi)=I(\psi;\eta)/I(\psi_0;\eta).
\]
Take rho0=central mass density, R=V0/sqrt(4pi G rho0), xi=r/R, and Q=-xi²psi'. Poisson equations:
\[
 \psi'=-Q/\xi^2,\quad Q'=\xi^2 x(\psi),\quad \psi(0)=\psi_0,\quad Q(0)=0. \tag{3}
\]
Let xi_t be first psi=0. Then finite total physical mass and circular velocity:
\[
 M_t=4\pi\rho_0 R^3 Q_t,\quad V_c^2=V_0^2 Q(\xi)/\xi. \tag{4}
\]
Potential and radial derivative match smoothly to vacuum Newtonian exterior, with additive potential constant fixing E_t.

## 3. Edge asymptotic derived analytically
For psi->0 at fixed eta>0, Taylor series of the BE occupation near E_t gives f≈exp(eta)/(exp(eta)-1)^2*(psi-q²). Hence
\[
 I(\psi;\eta)=\frac{2e^\eta}{15(e^\eta-1)^2}\,\psi^{5/2}[1+O(\psi)]. \tag{5}
\]
A finite Q_t makes psi near the edge linear in xi_t-xi, hence rho proportional to (xi_t-xi)^(5/2). This edge exponent also occurs in other lowered distributions, and is not claimed as new in statistical mechanics. Numerical exponents at eta=3,4,8 were 2.5000068, 2.5000064 and 2.5000061; leading coefficients agree at psi=1e-6 to 3e-7 relative.

## 4. Independent finite-radius solutions
At eta=4 with monotone interpolation of I and a high-accuracy potential-flux ODE integration:

| psi0 | xi_t=rt/R | Q_t=Mt/(4*pi*rho0*R³) | rho(0.5 rt)/rho0 |
|---:|---:|---:|---:|
| 1.0 | 5.9619843 | 2.1590945 | 0.054366 |
| 2.0 | 9.8649623 | 6.0747598 | 0.027279 |
| 3.0 | 15.6464877 | 11.4185729 | 0.009561 |
| 3.8 | 27.3724789 | 18.7511844 | 0.002216 |

For psi0=3.8 the central fugacity is z0=exp(psi0-eta)=exp(-0.2)≈0.819, NOT saturated z0=1. Therefore these solutions do not demonstrate a critically degenerate core from paper 1. The eta-psi0->0+ limit requires further analysis.

## 5. Tests, prior art and limitations
Ten validation checks passed: three positivity, three boundary asymptotic, four positive finite-radius and finite-mass calculations. These are numerical idealized equilibria, no galaxy velocity data, no cosmology, no stability or kinetic formation time. King/Woolley lowered isothermal distributions and fermionic King models are well-developed prior art, and the lowered prescription is not claimed original. New value relative to the *source model* is a concrete testable, positive truncated Bose variant with finite mass.

## 6. Next investigations
- Scan eta and psi0, include near-critical eta-psi0 approaching zero, then compare inner profiles to untruncated model.
- Evaluate radial Vlasov–Poisson perturbations and energy-Casimir sufficient stability criteria; do not assume all positive equilibria are stable.
- Fit the original 22 galaxies and independent baryon-aware SPARC samples, including the extra cutoff parameter penalty and held-out predictive accuracy.
- Derive a physical tidal/evaporative formation mechanism, a global cosmological and kinetic evolution and assess the effect of nonzero contact self interaction.

Sources: Abdullin & Popov (2021) https://doi.org/10.1088/1475-7516/2021/11/055; Chavanis review (2025) https://doi.org/10.3389/fspas.2025.1538434; classical King https://doi.org/10.1103/PhysRevD.91.063531; fermionic King https://doi.org/10.1103/PhysRevD.92.123527; tidal stripped halo DF https://doi.org/10.1093/mnras/stx652.