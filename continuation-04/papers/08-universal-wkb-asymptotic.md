# 08. Universal asymptotic *self-invalidation* of the local homogeneous Bose–Jeans scale near condensation

**Status:** analytically derived from original ideal Bose-Poisson equilibrium plus a marginal homogeneous response. No global stability theorem or observational exclusion is claimed.

## 1. Independent small parameter

Let \(\nu_c=\zeta(3/2)\), \(\Delta\nu=\nu_c-\nu_0>0\), \(\delta=-\log z_0>0\). As the central fugacity approaches saturation,
\[
\nu_0=\operatorname{Li}_{3/2}(e^{-\delta})
=\nu_c-2\sqrt\pi\sqrt\delta+O(\delta),
\qquad
\boxed{\delta=\frac{(\Delta\nu)^2}{4\pi}[1+O(\Delta\nu)]}.
\]
The central *collisionless homogeneous* zero-frequency susceptibility
\[
A(0)=\frac{\operatorname{Li}_{1/2}(e^{-\delta})}
{\operatorname{Li}_{3/2}(e^{-\delta})}
\sim\frac{\sqrt\pi}{\nu_c\sqrt\delta}
\sim\frac{2\pi}{\nu_c\Delta\nu}.
\]
Consequently the local homogeneous Jeans wavelength for \(g=0\) is
\[
\boxed{\lambda_J(0)=2\pi R A(0)^{-1/2}
\sim \sqrt{2\pi\nu_c}\,R\,(\Delta\nu)^{1/2}.}\tag{1}
\]
This result treats the central medium as infinite and uniform. It is an intermediate diagnostic, *not* a physical eigenmode of the galaxy.

## 2. Size of the region where that susceptibility exists

The critical effect only occurs near the galactic center. For the actual spherical stationary Poisson–Bose equilibrium and pure noncondensed component,
\[
y(r)=\frac{1}{6}(r/R)^2+O((r/R)^4)
\]
whenever \(\nu_0<\nu_c\), as \(r\to0\). Define the first half-response radius \(r_{1/2}\) by \(A(r_{1/2})=A(0)/2\), using \(z(r)=z_0e^{-y(r)}\). To leading order \(A(r)\propto(\delta+y)^{-1/2}\), so
\[
y(r_{1/2})\sim3\delta,\qquad
\boxed{r_{1/2}\sim R\sqrt{18\delta}
\sim \frac{3}{\sqrt{2\pi}}R\,\Delta\nu.}\tag{2}
\]
Thus **the high-response spatial region shrinks linearly in \(\Delta\nu\)**, whereas the homogeneous Jeans wavelength shrinks only as \(\sqrt{\Delta\nu}\). Dividing (1)/(2):
\[
\boxed{\frac{\lambda_J(0)}{r_{1/2}}
\sim \frac{2\pi\sqrt{\nu_c}}{3}\frac{1}{\sqrt{\Delta\nu}}
=3.3851381878\,(\Delta\nu)^{-1/2}
\longrightarrow\infty.}\tag{3}
\]

This is independent of halo scale \(R\), reference \(V_0\), chosen particle mass and the arbitrary comparison galaxy, *within the ideal static model*. A more precise definition of "response radius" modifies the order-one prefactor but not the divergent exponent. Adding repulsive HF interactions introduces a separate inner cutoff, analyzed in continuation-03 and note 06 of this continuation.

## 3. Quantitative agreement with direct polylogarithmic roots

The independent full special-function calculation uses the exact \(\nu(z)\), exact \(A(z)\), and the source's leading central potential, not the high-degeneracy expansion for the numerical column.

| \(\nu_0\) | Exact-root \(\lambda_J/r_{1/2}\) | Asymptotic (3) | asymptotic/exact |
|---:|---:|---:|---:|
| 2.2 | 4.49658 | 5.27145 | 1.1723 |
| 2.5 | 9.71336 | 10.09813 | 1.0396 |
| 2.6 | 30.30380 | 30.42970 | 1.00415 |
| 2.61 | 69.40139 | 69.45647 | 1.00079 |
| 2.612 | 174.70468 | 174.72657 | 1.00013 |

Convergence validates the \((\Delta\nu)^{-1/2}\) law. For \(\nu=2.2\), the central quadratic-potential approximation reaches large r/R and is rough; close-to-critical rows with \(r_{1/2}\ll R\) are the controlled limit.

## 4. What this rules out and what it emphatically does not

Equation (3) establishes a **mathematical failure of using the locally homogeneous Jeans plane-wave approximation ever more deeply into the saturated central cusp**. The candidate homogeneous marginal wavelength increasingly samples a strongly inhomogeneous background. It does not follow that the global halo is stable: an unstable non-local spherical eigenmode may exist, but it cannot be identified by attaching a universal local Jeans wavelength to the tiny high-A region.

This resolves the apparent paradox "critical Bose susceptibility diverges, therefore arbitrarily small clouds collapse": a local susceptibility can diverge even while the region with that susceptibility becomes vanishingly narrow. Only a global orbital/spectral treatment can establish the actual instability. The conclusion also limits naive identification of the HF ~377 AU local homogeneous scale with a Bose-star formation radius.

## 5. Next calculational step

Formulate global radial linearized Vlasov–Poisson response on finite-mass lowered Bose equilibrium, with finite-energy cutoff and an explicit perturbation class. Validate numerical spectral method on classical King models before varying \(\nu\) toward criticality and adding finite HF \(g\) or a coherent component.

**Prior art:** the general WKB requirement that wavelength be small against background variation is standard. The *particular asymptotic ratio (3)* is derived here for the Abdullin–Popov stationary ideal critical Bose model; novelty relative to the full kinetic literature has not been checked. Relevant sources: JCAP 2021 https://doi.org/10.1088/1475-7516/2021/11/055; 2023 Vlasov stability review https://doi.org/10.1088/1361-6382/acf436; 2022 Bose cluster equilibrium https://doi.org/10.1103/PhysRevE.105.034145.
