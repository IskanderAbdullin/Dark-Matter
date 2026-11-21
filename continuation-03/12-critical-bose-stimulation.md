# A critical singularity of stimulated occupation averages: finite moments, divergent susceptibility

**Date:** 2026-11-21. This is an analytic result of textbook Bose statistics applied explicitly to the JCAP 2021 halo distribution. Its algebraic identities are *not* claimed new; the physical use is as a rigorous warning against a uniform Bose-enhancement multiplier.

## 1. Exact identity for the ideal Bose distribution
Let \(f(p)=[z^{-1}\exp(p^2/2mkT)-1]^{-1}\), \(0<z<1\), and define \(\lambda_{\rm th}=\sqrt{2\pi}\hbar/(mV_0)\), \(V_0^2=kT/m\). In the usual quantum-cell measure \(d^3p/(2\pi\hbar)^3\),
\[
n=\lambda_{\rm th}^{-3}\,\mathrm{Li}_{3/2}(z).
\]
Expanding \(f=\sum_{s=1}^\infty z^s e^{-s\beta\epsilon}\) and squaring, the coefficient of the \(s\)-th term in \(f^2\) is \(s-1\), so
\[
\int\frac{d^3p}{(2\pi\hbar)^3}f(p)^2
=\lambda_{\rm th}^{-3}
\sum_{s=2}^\infty(s-1)\frac{z^s}{s^{3/2}}
=\lambda_{\rm th}^{-3}
\left[\mathrm{Li}_{1/2}(z)-\mathrm{Li}_{3/2}(z)\right].
\]
Therefore **the number-weighted mean multiplicative occupation factor** in the *same momentum cell* is
\[
\boxed{\langle1+f\rangle_n
=\frac{\int f(1+f)d^3p}{\int f\,d^3p}
=\frac{\mathrm{Li}_{1/2}(z)}{\mathrm{Li}_{3/2}(z)}}.\tag{1}
\]
The ratio is finite for every \(z<1\) and diverges as \(z\uparrow1\). This is a mathematically well-defined illustrative average, **not** the true collision-integrated enhancement of scattering; real \(2\leftrightarrow2\) rates weight initial and final momenta, energy-momentum conservation, matrix elements, inverse processes and occupation factors for different modes.

## 2. Critical asymptotic and connection to Hartree–Fock
Define \(\delta=-\ln z>0\), \(\nu=\mathrm{Li}_{3/2}(e^{-\delta})\) and \(\nu_c=\zeta(3/2)\). As \(\delta\to0^+\),
\[
\mathrm{Li}_{1/2}(e^{-\delta})
=\sqrt\pi\,\delta^{-1/2}+\zeta(1/2)+O(\delta^{1/2}),
\]
\[
\nu_c-\nu=2\sqrt\pi\,\delta^{1/2}+O(\delta).
\]
Combining with (1),
\[
\boxed{\langle1+f\rangle_n\sim
\frac{\sqrt\pi}{\zeta(3/2)\sqrt\delta}
\sim\frac{2\pi}{\zeta(3/2)\,[\nu_c-\nu]}}\tag{2}
\]
with coefficient \(2\pi/\zeta(3/2)\simeq2.40516\).
The **isothermal density susceptibility** is also proportional to the singular polylog:
\[
\left(\frac{\partial n}{\partial\mu_{\rm eff}}\right)_T
=\frac{\lambda_{\rm th}^{-3}}{kT}\mathrm{Li}_{1/2}(z).
\]
Thus the same susceptibility which makes an arbitrarily small repulsive Hartree–Fock shift \(2gn\) important at the critical mathematical center is the one making the weighted factor in (1) divergent. This analytically connects the earlier manuscript 6 finite-coupling boundary-layer result with the need for careful Bose-enhanced collisions.

## 3. Independent numerical checks
Exact mpmath high-precision polylogs and the asymptotic expansion were independently compared:

| \(\delta=-\ln z\) | \(\nu\) | \(\langle 1+f\rangle_n\) | \((\nu_c-\nu)\langle1+f\rangle_n\) |
|---:|---:|---:|---:|
| 0.1 | 1.6363774 | 2.54544 | 2.48434 |
| 0.01 | 2.2724777 | 7.15794 | 2.43297 |
| 0.0001 | 2.5770723 | 68.2111 | 2.40806 |
| 0.000001 | 2.6088319 | 678.8454 | 2.40545 |
| 0.00000001 | 2.6120209 | 6785.198 | 2.40519 |

Approaches \(2\pi/\nu_c=2.40516176603\) from above as expected.

## 4. Why this matters despite finite particle fractions
Our prior numerical estimate for an *exact critical ideal* gas: only ~26.25% of particles occupy modes with \(f>10\) and ~8.60% modes with \(f>100\). Nevertheless the **number-weighted mean** \(1+f\) diverges as the infrared cutoff is removed: a small tail of exceptionally occupied modes controls this particular moment.

This does **not** imply an infinite physical scattering rate. At exactly \(z=1\), \(f(p)\sim p^{-2}\) while \(p^2dp\,f(p)\) is integrable (finite n), but \(p^2dp\,f(p)^2\sim dp/p^2\) is not. Full collision integrals include energy-momentum constraints, gain-loss cancellations and possible finite-volume, condensate, interactions and dynamical effects that can regularize or restructure the IR behavior. Simply replacing \(\sigma/m\) by \((\sigma/m)\times(\text{huge occupancy})\) without deriving the process-specific kernel is unjustified.

## 5. Falsifiable immediate next computation
- Pick an explicit conserved-number contact Lagrangian and quantum Boltzmann \(C_{22}\) (or the Proukakis–Rigopoulos–Soto GP–Boltzmann system including \(C_{12}\)).
- Construct \(z_0\) values close to, but not exactly, unity, and compare collision rates at fixed local \(n,T\), enforcing conservation of mass, energy and detailed balance.
- Add a physically justified IR regulator (finite halo, self-coupling, condensate threshold), vary it, and test converged transport relaxation time rather than a formal occupancy moment.
- Couple resulting rates to the formation of a bound Bose-star population and to the cosmological initial distribution.
- Do not claim a specific observational or kinetic bound until the full kernel is evaluated.

**Provenance:** Abdullin & Popov 2021 (distribBE, nonconddens, ncEq); Proukakis–Rigopoulos–Soto PRD 108, 083513 (2023) and PRD 111, 023505 (2025), https://doi.org/10.1103/PhysRevD.108.083513, https://doi.org/10.1103/PhysRevD.111.023505; standard Bose polylog asymptotics. This note is a calculation and theoretical warning, not a claimed new finding in general statistical mechanics.
