# Phase-space density of states at the tidal boundary of a finite Bose halo

**21 November 2026 — independent analytic and numerical study.** Continues the critical finite-radius lowered-Bose equilibria in continuation-04 and the toy open-energy flux calculation in continuation-05. This is a STATIC phase-space measure and **does not** predict the physical quantum collision or evaporation coefficient.

## 1. Derive the orbit-energy state measure before specifying a Fokker–Planck escape flux

For an isotropic, spherical, time-independent equilibrium use the same dimensionless variables as the earlier Bose–Poisson analyses:

\[
\xi=r/R,\quad q=p/\sqrt{2mk_BT},\quad
\psi(\xi)=\frac{E_t-m\Phi(\xi)}{k_BT},\quad
s=\frac{E_t-E}{k_BT}=\psi-q^2.
\]

The particles of binding energy \(s\) occupy positions satisfying \(\psi(\xi)>s\), with outer turning radius \(\xi_s\) defined by \(\psi(\xi_s)=s\). Differentiating the six-dimensional phase-space volume with respect to **physical** single-particle energy gives

\[
\boxed{\frac{d\Gamma}{dE}=
\frac{8\pi^2 R^3(2mk_BT)^{3/2}}{k_BT}\,
\mathcal G(s),\quad
\mathcal G(s)=\int_0^{\xi_s}\xi^2\sqrt{\psi(\xi)-s}\,d\xi.}
\tag{1}
\]

The prefactor counts physical phase-space volume; divide it by \((2\pi\hbar)^3\) to count quantum cells, consistently with a dimensionless Bose occupation \(f(E)\). The energy distribution of particles is

\[
\frac{dN}{dE}=\frac{1}{(2\pi\hbar)^3}\frac{d\Gamma}{dE}\,f(E).
\tag{2}
\]

This quantity is not the spatial density \(\rho(r)\): two powers describing their edges need not be equal.

## 2. The escape state measure is finite and nonzero

For a positive self-consistent finite-radius halo with \(\psi(0)=\psi_0>0\), \(\psi(\xi_t)=0\), \(Q_t>0\), and \(\psi'(\xi_t)=-Q_t/\xi_t^2<0\), the edge potential is linear:

\[
\psi(\xi)=\frac{Q_t}{\xi_t^2}(\xi_t-\xi)+O((\xi_t-\xi)^2).
\]

Consequently \(0<\mathcal G(0)<\infty\), and the first derivative exists:

\[
\boxed{\mathcal G'(0)=-
\frac12\int_0^{\xi_t}
\frac{\xi^2\,d\xi}{\sqrt{\psi(\xi)}}<0.}
\tag{3}
\]

The integrable \((\xi_t-\xi)^{-1/2}\) endpoint is verified using a quadratic coordinate transformation in the numerical quadrature. Thus the energy-space state measure does **not** introduce a singular power at the escape boundary of these finite-radius examples.

If \(f(E)\sim C(E_t-E)^k\), then \(\frac{dN}{dE}\sim C\frac{d\Gamma}{dE}|_{E_t}(E_t-E)^k\), with the SAME exponent k, while the *spatial density* behaves \(\rho(r)\propto\psi(r)^{k+3/2}\). The 3/2 difference is the local velocity integral, not a divergence in the global energy-state measure.

## 3. Comparison of two previously verified finite Bose equilibria

The two models below have the SAME dimensionless central depth \(\eta=\psi_0=4\), \(R\), \(V_0\), and \(\rho_0\), but different assumed energy cutoffs. These are not fitted galaxies and do **not** have the same microscopic particle mass if full quantum normalization is imposed.

| Dimensionless observable | Linear lowering, k=1 | Quadratic C¹ lowering, k=2 |
|---|---:|---:|
| Outer radius \(\xi_t\) | 40.298792 | 130.178357 |
| Dimensionless mass \(Q_t\) | 26.357019 | 26.319405 |
| Edge energy density of states \(\mathcal G(0)\) | 10,194.818 | 194,152.593 |
| \(\mathcal G'(0)\) | −39,794.736 | −2,408,749.823 |
| \(\mathcal G'(0)/\mathcal G(0)\) | −3.90343 | −12.40648 |

The quadratic-lowered model has only ~0.143% lower mass, yet **19.044 times more one-particle state-volume at the escape energy**. This is a strong conditional demonstration that the dilute outer geometry can radically alter the near-escape phase-space capacity without substantially changing total mass.

**Important:** a factor of 19 in \(\mathcal G(0)\) is **not** a factor of 19 in the true escape rate. Orbit-averaged diffusion, potential evolution, particle density, Bose stimulation, mass flux, and tidal conditions enter the full kinetic operator separately. This calculation only isolates one geometry factor.

Direct numerical integrals at \(s=10^{-4}\) and \(10^{-3}\) reproduce the independently integrated \(\mathcal G'(0)\) within 1% for both profiles. Positive monotonic \(\mathcal G(s)\) is observed at the sampled energies.

## 4. Universal shallow-cutoff theorem

The previous exact-critical lowered-Bose family admits a controlled \(\eta\to0^+\) limit (continuation-04 paper 11):

\[
\psi(\xi)=\eta U(u)+o(\eta),\quad \xi=\sqrt\eta\,u,\quad
u_t=11.860598,\quad \mathcal Q_t=2.895576.
\]

At *fixed scaled binding energy* \(s=\eta\,z\), substitute this rescaling in (1). Each spatial factor gives \(\xi^2d\xi=\eta^{3/2}u^2du\), while the velocity shell gives \(\sqrt{\psi-s}=\sqrt\eta\sqrt{U-z}\). Hence

\[
\boxed{\mathcal G_\eta(\eta z)=\eta^2\,\mathfrak G(z)+o(\eta^2),\qquad
\mathfrak G(z)=\int_0^{u_z}u^2\sqrt{U(u)-z}\,du.}
\tag{4}
\]

In particular, direct solution of the universal Poisson ODE and an independent energy-shell integral yield

\[
\boxed{\mathfrak G(0)=157.329102882.}
\tag{5}
\]

An independently integrated finite-\(\eta\) sequence approaches this constant:

| \(\eta\) | \(\mathcal G_\eta(0)/\eta^2\) | Relative deviation |
|---:|---:|---:|
| 0.5 | 160.656034 | 2.115% |
| 0.2 | 157.856700 | 0.335% |
| 0.1 | 157.461005 | 0.0838% |
| 0.05 | 157.362354 | 0.0211% |
| 0.02 | 157.334557 | 0.00347% |
| 0.01 | 157.330598 | 0.000950% |

Equation (4) is a statement about **static equilibria under the chosen nondimensional normalization**, not about the rate at which a real self-gravitating gas loses particles.

## 5. Implications for the next physical program

1. Implement the correctly normalized energy-space particle number \(N(E)=g(E)f(E)\) in an orbit-averaged kinetic solver; check number and energy conservation with an explicit escape boundary.
2. Derive the energy-dependent diffusion and friction *from the actual Bose \(2\leftrightarrow2\) collision kernel and gravitational encounters* rather than setting a constant coefficient.
3. Check whether a linear edge \(f(E)\propto(E_t-E)\), as selected by the simple finite-D open-flux toy, persists after including the state density and real escape geometry.
4. Distinguish the **state-counting** function in (1) from a dynamical relaxation time and the cosmological linear transfer function; none is computed here.
5. Refit realistic tidal radii and stellar nuisance parameters before interpreting the outer geometric factor as an observational signal.

### Literature and prior art

The use of energy-space density of states in orbit-averaged Fokker–Planck models is established: [orbit-averaged energy-space transport review](https://doi.org/10.1007/s41114-018-0013-8) defines \(N(E)=16\pi^2p(E)F(E)\). [Self-similar orbit-averaged Fokker–Planck models](https://doi.org/10.1088/1674-4527/21/5/108) explicitly note the linear King distribution near escape. A new 2026 mass-preserving numerical construction of Vlasov steady states offers a separate future benchmark, https://doi.org/10.1016/j.nonrwa.2025.104467. The **formulas for this specific Bose equilibrium and the numerical geometric comparison** are calculations in this research project, not claims of invention of the energy-density-of-states method.
