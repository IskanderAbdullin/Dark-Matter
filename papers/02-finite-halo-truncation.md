# Paper 2 — A finite-mass lowered Bose–Einstein halo with compact condensate populations
**Research proposal and analytic model, 2026-11-21.** No numerical stability or observational fit is claimed.
**Working title:** *Energy-truncated noncondensed bosonic dark-matter equilibria and self-consistent condensed substructure*

## Abstract
The noncondensed equilibrium Bose halo of Abdullin & Popov (2021) approaches an isothermal-sphere envelope at low occupation. Such an isolated isothermal envelope has asymptotically divergent total mass. We propose a positive, energy-truncated or “lowered” Bose distribution, derive the equations needed to obtain finite-radius equilibria, and couple its gravitational potential to a mass-weighted population of bound Bose-star substructures. The goal is a predictive extension whose virial boundary and total mass are explicit, enabling lensing, cosmological and galaxy-to-galaxy comparisons that cannot be made with an unbounded isothermal solution.

## 1. Why the untruncated halo is insufficient
For small fugacity the Bose distribution tends to Maxwell–Boltzmann statistics. A singular isothermal asymptote has \(\rho(r)\sim V_0^2/(2\pi G r^2)\). Therefore \(M(<r)\propto r\) and its asymptotically flat rotation curve has \(V_\infty\sim\sqrt2 V_0\). A finite virial radius drawn on this solution is an external cutoff, not a stationary equilibrium boundary. A proper comparison to cosmological mass–concentration relations or wide-field lensing requires a physical prescription at escape energy.

## 2. Candidate lowered Bose distribution
Write \(E=p^2/(2m)+m\Phi(r)\), and let \(E_t\) be the bound/unbound escape-energy threshold. For \(\mu<m\Phi(0)\) define
\[
F(E)=\frac{\Theta(E_t-E)}{(2\pi\hbar)^3}
\left\{\frac{1}{e^{\beta(E-\mu)}-1}
-\frac{1}{e^{\beta(E_t-\mu)}-1}\right\},
\qquad \beta=(k_BT)^{-1}. \tag{1}
\]
This distribution is nonnegative for \(E<E_t\), vanishes continuously at the boundary, and approaches the ordinary Bose function far below the escape threshold if the subtraction term is negligible. This is **a proposed collisionally relaxed distribution**, not necessarily the outcome of collisionless violent relaxation.

At radius \(r\) let \(\Psi(r)=E_t-m\Phi(r)>0\) be the relative escape energy. The local momentum cutoff is
\[
0\le p\le \sqrt{2m\Psi(r)}.
\]
The number density and scalar pressure become finite integrals:
\[
n_{\rm n}(r)=4\pi\int_0^{\sqrt{2m\Psi}}F(E)\,p^2\,dp,\quad
P_{\rm n}(r)=\frac{4\pi}{3m}\int_0^{\sqrt{2m\Psi}}F(E)\,p^4\,dp. \tag{2}
\]
Dimensionless quadrature of (2), not evaluation of the untruncated polylog, must be used at each potential step.

## 3. Finite-mass boundary asymptotic
Let \(\Delta=E_t-E\) be small. Because \(F(E_t)=0\), a regular Taylor expansion gives \(F(E)\simeq C(E_t-E)\), where
\[
C=\frac{\beta e^{\beta(E_t-\mu)}}{(2\pi\hbar)^3[e^{\beta(E_t-\mu)}-1]^2}>0.
\]
Near a boundary \(\Psi\downarrow0\), integration over the allowed momenta gives
\[
n_{\rm n}\propto\Psi^{5/2},\qquad
P_{\rm n}\propto\Psi^{7/2}. \tag{3}
\]
This establishes a finite density at the edge and a vanishing boundary pressure, assuming a finite tidal radius exists. It does not by itself demonstrate a stable equilibrium or unique total halo mass.

## 4. Self-consistent Poisson and substructure
For the smooth boson component plus a mean condensed-object component,
\[
\frac1{r^2}\frac{d}{dr}\!\left(r^2\frac{d\Phi}{dr}\right)
=4\pi G\left[m n_{\rm n}(\Psi)+\rho_{\rm BS}(r)+\rho_{\rm b}(r)\right]. \tag{4}
\]
The condensate population is a *mass function*, not an arbitrary extra density:
\[
\rho_{\rm BS}(r)=\int_{M_{\min}}^{M_{\max}}M\,
\frac{d^2N_{\rm BS}}{dV\,dM}(r,M)\,dM. \tag{5}
\]
Its isotropic Jeans limit, when admissible, gives \(\rho_{\rm BS}(r)\propto \exp[-(\Phi(r)-\Phi_0)/v_{\rm BS}^2]\); for finite mass this must likewise be truncated at an escape surface. The two central parameters \(\alpha=\rho_{\rm BS}(0)/\rho_{\rm n}(0)\) and \(\gamma=V_0^2/v_{\rm BS}^2\) are independent physical variables unless a dynamical mechanism relates them. This explicitly relaxes the source paper's simplifying \(\alpha\simeq\gamma\) choice.

## 5. Boundary-value problem and numerical programme
Unknowns: the potential and its derivative, fugacity/chemical potential, velocity scale \(V_0\), escape threshold \(E_t\), condensed fraction, stellar/gas baryonic fields, and external tidal/virial conditions. Use central regularity \(\Phi'(0)=0\), \(0<z_0<1\), integrate Poisson while updating \(n,P\) from (2), then require \(\Psi(r_t)=0\). There may be multiple branches. Solve by shooting/continuation, validate convergence near \(z_0\to1\), and explore Antonov/gravothermal instability by a proper entropy and energy analysis.

Suggested numerical tests: (i) positivity \(F\ge0\); (ii) mass closure \(M'(r)=4\pi r^2\rho\); (iii) dependence on \(E_t\); (iv) conventional lowered Maxwellian limit for \(z_0\to0\); (v) recovery of the untruncated interior as \(\Psi(0)\) grows; (vi) robustness of lensing convergence \(\Sigma(R)\). A single stable solution is not enough to establish a population prior.

## 6. Discriminating predictions
The lowered model predicts a profile break in the outskirts linked to escape energy and tidal environment. Joint rotation curves plus outer satellite kinematics or weak lensing might constrain that break. Unlike ad hoc cutoff, the lowered model enforces phase-space positivity and zero edge pressure. If no stable boundary solution can fit observed outer-profile behavior with a common boson mass and acceptable baryonic content, the extension fails.

## 7. Open risks and citations
The distribution is a modeling ansatz; its thermodynamic derivation and possible distinction from a violent-relaxation distribution require work. Baryonic response, anisotropic orbit populations, tides, and time-dependent accretion must eventually be included. Relevant precedents: [Abdullin & Popov (2021)](https://doi.org/10.1088/1475-7516/2021/11/055); [Chavanis (2025) review](https://doi.org/10.3389/fspas.2025.1538434); Lynden-Bell (1967) on violent relaxation; standard lowered-isothermal King models (must be separately surveyed and cited before publication).


## Appendix — numerical verification of the infinite-mass isothermal envelope
A [new independent log-radius integration](../continuation-02/10-collisionality-and-asymptotic-tail.md) verifies the outer asymptote of the *untruncated* critical Bose halo, \(\rho/\rho_0\sim2(R/r)^2\), \(u^2=V^2/V_0^2\to2\), giving \(M(<r)\propto r\) without a finite total mass. Numerically, at \(\xi=r/R=10^5\), \(\xi^2 x=2.00508\) and \(u^2=2.00053\). These are static equilibrium results at very large mathematical radii, not a claim that an actual halo extends to \(10^5R\). A physical tidal/escape-energy truncation must be implemented consistently in the distribution function and Poisson potential, rather than merely multiplying the mass density by an arbitrary radial window.


## Appendix B — solved positive finite-mass lowered Bose distribution
The [new calculation](../continuation-02/11-lowered-bose-finite-halo.md) uses a *phase-space* escape-energy cutoff rather than a radial density taper. For energy E below E_t, select
\[
 f(E)=\frac{1}{e^{(E-\mu)/k_BT}-1}-\frac{1}{e^{(E_t-\mu)/k_BT}-1};
 \quad f=0 \text{ for } E\ge E_t .
\]
This is nonnegative when E>mu and a stationary collisionless Vlasov–Poisson ansatz. Its edge density scales as rho~psi^(5/2), where psi=(E_t-mPhi)/(k_BT). At eta=(E_t-mu)/(k_BT)=4, numerical solutions have finite radii rt/R=5.962,9.865,15.646,27.372 and dimensionless masses Q_t=2.159,6.075,11.419,18.751 for psi0=1,2,3,3.8. Ten positivity/asymptotic/existence checks passed. This is not a proof of thermodynamic equilibrium, dynamical stability or an observed-galaxy fit. Requiring psi0<eta forces z0=exp(psi0-eta)<1, so this demonstration does not include the exactly saturated central fugacity.
