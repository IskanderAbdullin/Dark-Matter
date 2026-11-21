# 07. Exact pressure–density derivative and independent virial verification of the **critical finite-mass lowered-Bose halo**

**Date:** 2026-11-21. This deepens paper 04 of the current continuation. The result is a self-consistency proof and high-accuracy numerical test *within the stationary spherical energy-cutoff ansatz*, **not** a stability theorem, a thermalization calculation, or a fit to observations.

## 1. Why the virial theorem is a useful additional independent test

The proposed lowered Bose occupation is \(f(q;\psi,\eta)= [e^{\eta-\psi+q^2}-1]^{-1}-[e^{\eta}-1]^{-1}\) for \(0\le q\le\sqrt{\psi}\), zero otherwise. With \(\psi=(E_t-m\Phi)/(kT)\), \(\eta=(E_t-\mu)/(kT)\), \(q=p/\sqrt{2mkT}\), define the density moment
\[
I(\psi,\eta)=\int_0^{\sqrt\psi}q^2 f(q;\psi,\eta)\,dq
\]
and the pressure moment
\[
J(\psi,\eta)=\int_0^{\sqrt\psi}q^4 f(q;\psi,\eta)\,dq .
\]
At the source's central density normalization \(\rho_0\), \(I_c=I(\psi_0,\eta)\), \(R^2=V_0^2/(4\pi G\rho_0)\), the **actual isotropic velocity dispersion** and pressure are
\[
\rho=\rho_0 I/I_c,\qquad
P=\frac23\rho_0V_0^2\,J/I_c,\qquad
\sigma_{1D}^2=\frac{2}{3}V_0^2 J/I.
\]
These differ from a globally fixed \(P=\rho V_0^2\) when degenerate or energy-lowered.

## 2. An exact identity rather than a postulated polytropic equation of state

At the upper momentum limit \(q=\sqrt\psi\), the subtracted occupation is **zero**. Its derivative with respect to \(\psi\) is \(\partial_\psi f=f_{\rm BE}(1+f_{\rm BE})\), whereas \(d_q f=-2q\,\partial_\psi f\). Integrating by parts and using vanishing boundary occupation,
\[
\frac{dJ}{d\psi}=\int_0^{\sqrt\psi}q^4\frac{\partial f}{\partial\psi}dq
=-\frac12\int_0^{\sqrt\psi}q^3\frac{\partial f}{\partial q}dq
=\frac32 I(\psi).
\]
Consequently
\[
\boxed{\frac{d}{d\psi}\left[\frac23 J(\psi)\right]=I(\psi)}.\tag{1}
\]
Since \(d\psi/dr=-m\Phi'/kT\), (1) gives exactly \(dP/dr=-\rho\Phi'\): the kinetic isotropic Jeans / hydrostatic balance follows from the same positive distribution used in Poisson. At the escape boundary \(I,J\to0\); no artificial finite boundary pressure is needed.

## 3. Scalar virial identity for a finite-radius equilibrium

Define \(x=I/I_c\), \(p=(2/3)J/I_c\), and \(Q(\xi)=G M(<r)/(V_0^2 R)\), \(\xi=r/R\). For a stationary, non-streaming, isolated, spherically symmetric self-gravitating halo,
\[
K=\frac32\int P\,d^3r,\qquad W=-4\pi G\int_0^{r_t}\rho M r\,dr,
\]
so in units \(C\equiv4\pi\rho_0V_0^2R^3\)
\[
\frac{K}{C}=\frac32\int_0^{\xi_t}p(\xi)\xi^2\,d\xi
=\int_0^{\xi_t}\frac{J}{I_c}\xi^2d\xi,\quad
\frac{W}{C}=-\int_0^{\xi_t}x(\xi)Q(\xi)\xi\,d\xi.
\]
Using \(\mathrm dP/\mathrm dr=-\rho\Phi'\) and \(P(r_t)=0\), one gets \(\boxed{2K+W=0}\). Unlike an arbitrary density cutoff, the chosen phase-space cutoff satisfies this identity as a consequence of the stationary Jeans theorem and Poisson equation.

## 4. Numerically verified critical endpoint \(\eta=\psi_0=4\)

Independent Gauss–Legendre quadrature for \(J\), reusing the Poisson solution and an independent radial integral, gave

| Quantity | Value |
|---|---:|
| Edge \(\xi_t=r_t/R\) | 40.2987920 |
| Mass \(Q_t\) | 26.3570198 |
| \(K/C\) | 17.503121767 |
| \(W/C\) | −35.006243536 |
| \((2K+W)/|W|\) | **−4.13×10⁻¹¹** |
| Maximum relative numerical error in \(d[(2/3)J/I_c]/d\psi=I/I_c\) at 7 interior probe radii | **1.44×10⁻⁶** |
| Normalized edge pressure | ~10⁻⁶¹ |

This is a useful high-precision **stationary-DF sanity check**. It does not prove that a Maxwellian-like or BE phase-space state is dynamically selected from cosmological initial conditions. A negative potential term in the second variation of the energy–Casimir functional is still present, so virialization and monotone \(df/dE\) are insufficient for a rigorous global nonlinear stability claim.

## 5. Further theoretical consequences

1. The new model has a *radially changing* local actual dispersion \(\sigma_{1D}^2=(2/3)V_0^2 J/I\) even though the occupation is parameterized by fixed \(T\). Comparing with the original untruncated \(\sigma_{1D}^2=V_0^2\mathrm{Li}_{5/2}/\mathrm{Li}_{3/2}\) gives an observable-shape and Jeans-consistency diagnostic.
2. The \(J\)-derivative identity provides a stronger code test for future cutoff/interactions implementations than checking only positivity and finite mass.
3. To assess gravothermal stability, compute entropy/energy curves, turning points and perturbation spectra with the *correct* lowered Bose quantum kinetic entropy; equilibrium virial is necessary but not sufficient.
4. To attach realistic halos, replace the vacuum exterior by a tidal environment and check escaping populations; a finite-radius standalone solution need not be globally cosmologically realized.

Relevant prior work: Abdullin & Popov (2021) https://doi.org/10.1088/1475-7516/2021/11/055; Kirejczyk et al. (2022) https://doi.org/10.1103/PhysRevE.105.034145; Guo & Rein (2007) https://doi.org/10.1007/s00220-007-0212-8; finite distribution truncations Baes (2026) https://doi.org/10.1093/mnras/stag050 and earlier https://doi.org/10.1093/mnras/stad117; Chavanis (2025) https://doi.org/10.3389/fspas.2025.1538434. The virial identity is **classical/known**, not claimed new in kinetic theory; the specific finite-mass critical numerical verification is new to this project.
