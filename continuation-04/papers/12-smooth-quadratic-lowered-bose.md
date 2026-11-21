# 12. A smooth, positive, finite-mass *quadratically lowered* Bose distribution

**New research proposal and numerical solution, 21 November 2026.** Prompted by the mathematical distinction between oscillating and damped Vlasov steady states with different vacuum-boundary regularity. No stability or Landau-damping theorem is claimed for this particular 3D isotropic Bose halo.

## 1. Motivation from direct stability theory (not a transfer of its theorem)

Hadžić, Rein, Schrecker & Straub, *Damping Versus Oscillations for a Gravitational Vlasov–Poisson System*, ARMA 249, 45 (2025), https://doi.org/10.1007/s00205-025-02114-y, prove a sharp regularity-dependent dichotomy in a **specific radial model with a fixed central point mass and fixed angular momentum**: for energy-cutoff polytropes \(f\propto(E_0-E)^k_+\), linear damping when \(k>1\) and undamped modes when \(1/2<k\le1\). Those hypotheses are **not satisfied by our 3D fully isotropic self-gravitating Bose halo**; hence their result is only a motivation to investigate smoothness.

Our previously constructed *linearly lowered* Bose distribution (continuation-02/04) satisfies \(f\sim{\rm const}\,(E_t-E)\) at the escape threshold, so it is continuous but not differentiable there (analog of \(k=1\)). We now define a mathematically positive C¹ variant with edge exponent \(k=2\), without claiming that the 2025 theorem guarantees damping.

## 2. Convex-tangent subtraction preserves nonnegative phase-space occupation

Write \(B(x)=(e^x-1)^{-1}\), \(x=(E-\mu)/(kT)>0\), and \(\eta=(E_t-\mu)/(kT)>0\). Let \(B'_t=B'(\eta)<0\) and \(c=B(\eta)\). For energies below escape \(E<E_t\), choose
\[
\boxed{
f_2(E)=B(x)-B(\eta)-B'(\eta)(x-\eta),
\qquad f_2(E)=0\quad (E\ge E_t).}\tag{1}
\]
For all \(x>0\), \(B''(x)=B(x)[1+B(x)][1+2B(x)]>0\). A convex function lies above every tangent, so
\[
\boxed{f_2(E)\ge0\quad\text{throughout its physical domain}.}
\]
At cutoff \(x\to\eta^-\), its Taylor expansion is
\[
f_2(E)=\tfrac12 B''(\eta)(\eta-x)^2+O((\eta-x)^3),
\]
hence \(f_2(E_t)=f_2'(E_t)=0\): unlike ordinary King-style linear lowering, the **piecewise** distribution is C¹ across the escape boundary. Its second derivative has a jump, so it is not C² there. For \(\eta=\psi_0\), the central ideal Bose singularity at \(p=0,r=0\) remains; the edge C¹ property does not make the entire distribution globally smooth.

In dimensionless \(q=p/\sqrt{2mkT}\), \(\psi=(E_t-m\Phi)/(kT)\), \(\Delta=\psi-q^2\), the formula is
\[
f_2(q;\psi,\eta)=f_1(q;\psi,\eta)-a_\eta(\psi-q^2),\quad
a_\eta=B(\eta)[1+B(\eta)]>0 .
\]

## 3. Exact density and pressure moments; stronger edge regularity

If \(I_1=\int_0^{\sqrt\psi}q^2f_1dq\) and \(J_1=\int_0^{\sqrt\psi}q^4f_1dq\) denote previous linear-lowered Bose moments, the new moments satisfy the **exact** identities
\[
\boxed{I_2(\psi;\eta)=I_1(\psi;\eta)-\frac{2a_\eta}{15}\psi^{5/2},}\tag{2}
\]
\[
\boxed{J_2(\psi;\eta)=J_1(\psi;\eta)-\frac{2a_\eta}{35}\psi^{7/2}.}\tag{3}
\]
As \(\psi\downarrow0\),
\[
\boxed{I_2(\psi;\eta)\sim\frac4{105}B''(\eta)\psi^{7/2}.}\tag{4}
\]
The edge density therefore scales as \(\rho\propto(r_t-r)^{7/2}\), instead of the \(5/2\) law for linearly lowered Bose. Both models have finite mass and vanishing edge pressure.

Since \(f_2\) depends only on the conserved particle energy and vanishes at escape, integration by parts gives the same exact isotropic Jeans identity \(dJ_2/d\psi=(3/2)I_2\). Poisson plus zero boundary pressure then implies \(2K+W=0\) in the stationary isolated model. This is a **consistency identity**, not proof of stability.

At exact central saturation \(\psi_0=\eta\), the singular \(-(\pi/2)\sqrt y\) reduction of the density moment remains unchanged in absolute normalization, but its normalized central radial coefficient is \(-x'(0^+)=\pi/(2\sqrt6 I_{2,c})\). Energy smoothing thus changes both outer falloff and central normalization.

## 4. Independently solved \(\eta=4\) equilibrium and virial

A new high-accuracy Poisson solver using the positive density (2) and carefully evaluated edge cancellation gives, for \(\eta=\psi_0=4\):

| Observable | Linear lowering \(f_1\) | C¹ quadratic lowering \(f_2\) |
|---|---:|---:|
| Central integral \(I_c\) | 1.0872575 | **1.00616755** |
| Phase normalization \(n_0\lambda^3\) | 2.45367744 | **2.27067701** |
| Central linear density slope \(-x'(0^+)\) | 0.589810 | **0.637344** |
| Outer density \(\rho\) | \(\propto\psi^{5/2}\) | \(\propto\psi^{7/2}\) |
| Cutoff radius \(r_t/R\) | 40.29879 | **130.17836** |
| Dimensionless total mass \(Q_t\) | 26.35702 | **26.31941** |

The radius changes by a factor \(\simeq3.23\) at fixed *central* \(\rho_0,V_0,\eta\), while total dimensionless mass remains quite similar. This large change in the dilute outer radial support emphasizes that finite-radius predictions depend materially on *how* the occupation is lowered; it is not a measurement of an actual galaxy size.

A direct momentum-integral edge exponent fit returned **3.500875** (target \(7/2\)). Direct positivity and central-slope assertions passed. With \(J_2\) calculated independently through (3), a separate high-resolution spherical integral returned
\[
\boxed{K/(4\pi\rho_0V_0^2R^3)=14.36038236,\qquad
W/(4\pi\rho_0V_0^2R^3)=-28.72076473,}
\]
so \((2K+W)/|W|\simeq+3.17\times10^{-11}\). This validates the stationary distribution and numerical moments, **not nonlinear perturbation stability**.

At the *same arbitrarily transplanted* F583-4 values \(R=0.6\,\mathrm{kpc},V_0=44.1\,\mathrm{km/s}\), the correct mass normalization would correspond to \(m\simeq89.28\) eV rather than the linearly lowered ~87.57 eV. Neither is a fresh observed mass inference; the new halo shapes must be refitted to real galaxy data.

## 5. Research payoff and falsifying tests

1. Compare linear- and quadratic-lowered families at fixed observable inner mass and realistic tidal environment, not just at fixed internal \(R,V_0\).
2. Repeat the global linearized Vlasov–Poisson \(l=0,1,2\) calculation for both edge regularities; only then test whether the damped/oscillatory dichotomy suggested by the 2025 *different* central-mass/fixed-angular-momentum model has an analogue.
3. Quantify actual escape/evaporation kinetics: a smooth cut-off is a stationary *ansatz*, not a self-consistently derived time-dependent boundary.
4. Reassess the all-particle-retained thermal-relic mixing test (paper 10). Any **exact** sharp finite-support DF, including the new C¹ one, still fails the closed thermal-patch *low-F* support necessary condition; smoothing the edge's first derivative is insufficient to permit pure closed collisionless formation.
5. Evaluate thermodynamic and dynamical stability separately, with the proper quantum entropy and energy–Casimir perturbation class; positive \(f_2\), C¹ edge and virial balance do **not** prove stability.

**References:** Hadžić et al. 2025 https://doi.org/10.1007/s00205-025-02114-y; King-like lowered phase-space distributions and stability Guo & Rein 2007 https://doi.org/10.1007/s00220-007-0212-8; Abdullin & Popov JCAP 2021 055 https://doi.org/10.1088/1475-7516/2021/11/055.
