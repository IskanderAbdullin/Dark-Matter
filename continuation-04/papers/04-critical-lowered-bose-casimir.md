# 04. Exactly critical *lowered* Bose halo: finite radius and an explicit convex energy–Casimir

**Date:** 2026-11-21. This is a new quantitative continuation of manuscript 2/11, which had solved only the strictly subcritical \(0<\psi_0<\eta\) equilibria. It does not establish nonlinear kinetic stability or microphysical production.

## A. Physically admissible cutoff ansatz
Write the single-particle energy \(E=p^2/2m+m\Phi(r)\) and cutoff \(E_t\). With \(\eta=(E_t-\mu)/kT>0\), \(\psi=(E_t-m\Phi)/kT\), select \(0\le\psi\le\eta\) and
\[
f(E)=\begin{cases}
\left[e^{(E-\mu)/kT}-1\right]^{-1}
-\left[e^{(E_t-\mu)/kT}-1\right]^{-1},
&E<E_t,\\
0,&E\ge E_t .
\end{cases}
\]
It is nonnegative for \(E\in(\mu,E_t)\) and depends only on E, hence is a stationary collisionless Vlasov ansatz for a self-consistent Poisson potential. This is a Bose analogue of **existing** King-like lowered energy distributions, not an original general technique.

Set \(q=p/\sqrt{2mkT}\) and
\[
I(\psi;\eta)=\int_0^{\sqrt\psi}dq\,q^2
\left[\frac1{e^{\eta-\psi+q^2}-1}-\frac1{e^\eta-1}\right].
\]
Central density normalization \(x=I(\psi;\eta)/I(\psi_0;\eta)\), \(R^2=V_0^2/(4\pi G\rho_0)\), \(\xi=r/R\), \(Q=-\xi^2\psi'\) obey
\[
\psi'=-Q/\xi^2,\quad Q'=\xi^2x(\psi),\quad \psi(0)=\psi_0,\quad Q(0)=0.
\]
The first zero \(\xi_t\) of \(\psi\) is an escape boundary; \(M=4\pi\rho_0R^3Q_t\) is finite. The 2026-11-21 previous work solved \(\eta=4\), \(\psi_0=1,2,3,3.8\).

## B. A mathematically new endpoint **within our program**: exactly saturated central fugacity

Unlike subcritical \(\psi_0<\eta\), now set \(\psi_0=\eta\), i.e. \(z_0=1\) and \(f(p=0,r=0)\) diverges as an ideal occupation, but the **three-dimensional density integral remains finite**:
\[
I_c(\eta)=I(\eta;\eta)
=\int_0^{\sqrt\eta}\!dq\,q^2\left[\frac1{e^{q^2}-1}-\frac1{e^\eta-1}\right].
\]
For \(\psi=\eta-y\) and small \(y>0\), an integrable near-zero momentum expansion yields
\[
I(\eta-y;\eta)=I_c(\eta)-\frac{\pi}{2}\sqrt{y}+O(y).
\]
Near a regular gravitational center, \(y=\xi^2/6+O(\xi^3)\), so the critical lowered-Bose *spatial* density is
\[
\boxed{x(\xi)=1-\frac{\pi}{2\sqrt6 I_c(\eta)}\,\xi+O(\xi^2).}\tag{1}
\]
It has finite central density and a one-sided linear radial slope (not a singular \(1/r\) cusp). The pure untruncated critical coefficient in manuscript 1 is \(0.5539795458\); lowering changes the coefficient because \(I_c\) is smaller than the full Bose momentum integral.

At \(\eta=4\), \(I_c=1.0872575050\), and the slope in (1) is **0.5898096009**. Three radii \(0.002,0.001,0.0005R\) numerically returned slopes \(0.589452,0.589627,0.589762\), converging toward the analytic value. A separate direct integral reproduced \(\pi/2\) as the coefficient of \(\sqrt y\) at small y.

## C. First full critical finite-mass equilibrium in this project

Direct 512-node Gauss–Legendre quadrature, endpoint-clustered monotone interpolation and high-accuracy Poisson ODE with a zero-potential event yielded:
\[
\boxed{\eta=\psi_0=4:\quad r_t/R=40.2987921,\qquad
Q_t=M/(4\pi\rho_0R^3)=26.3570194.}\tag{2}
\]
The result is stable to changing quadrature nodes 96,192,256,512: \(r_t/R\) agrees to ~\(5\times10^{-6}\), \(Q_t\) to ~\(7\times10^{-7}\). No unphysical negative phase-space occupancy appears away from the central zero-momentum ideal singular limit.

**For illustration only**, transplanting \((R,V_0)=(0.6\,\mathrm{kpc},44.1\,\mathrm{km/s})\) from the **different untruncated model's** F583-4 fit would imply \(r_t\simeq24.18\) kpc and \(M_t\simeq7.15\times10^9\,M_\odot\). The correctly normalized phase density in this truncated profile is \(\nu_{\rm lowered}=4I_c/\sqrt\pi=2.4536774\); the mass corresponding to these *same assumed* R and V0 is \(m\simeq87.5663\) eV (not literally the original 86 eV). **These are not observational fits or credible virial estimates for F583-4.**

## D. Exact inverted distribution and a convex energy–Casimir density

For \(0<\psi_0<\eta\), let \(c=[e^\eta-1]^{-1}>0\). If \(f(E)= [e^{(E-\mu)/kT}-1]^{-1}-c\), then for \(f\ge0\)
\[
\frac{E(f)-\mu}{kT}=\log\left(1+\frac1{f+c}\right),\qquad
\frac{E_t-\mu}{kT}=\eta.
\]
Define \(\mathcal C'(f)=(E_t-E(f))/kT=\eta-\log(1+1/(f+c))\). The explicit convex primitive \(\mathcal C(0)=0\) is
\[
\mathcal C(f)=\eta f-\left[
(1+f+c)\log(1+f+c)-(f+c)\log(f+c)
-(1+c)\log(1+c)+c\log c
\right],
\]
with
\[
\boxed{\mathcal C''(f)=\frac{1}{(f+c)(1+f+c)}>0.}\tag{3}
\]
This makes the chosen \(f(E)\) a stationary point of the **formal** mass-constrained energy–Casimir functional
\[
\mathcal F[f]=H_{\rm VP}[f]+kT\!\int\! \mathcal C(f)\,d\Gamma
-E_t\!\int\! f\,d\Gamma.
\]
One must use a consistent phase-space mass/occupation normalization; in this formula H is particle-energy integral plus self-gravitating potential energy.

**Crucially, \(\mathcal C''>0\) does NOT prove \( \delta^2\mathcal F>0\)**, since gravity gives a negative, nonlocal contribution:
\[
\delta^2\mathcal F \sim \frac{kT}{2}\int\mathcal C''(f_0)(\delta f)^2 d\Gamma
-\frac{G}{2}\int d^3x\,d^3y\frac{\delta\rho(x)\delta\rho(y)}{|x-y|}.
\]
Actual theorems additionally require functional compactness, accessible perturbations, symmetry and possible cutoff regularity. At exactly critical center \(\sup f_0=\infty\), \(\mathcal C''(f_0)\to0\) pointwise, making a uniform-coercivity proof especially delicate. For an illustrative subcritical \(\eta=4,\psi_0=3.8\), \(c=0.01865736\), \(f_{\max}=4.497998\), \(\mathcal C''(f_{\max})\approx0.040134\); second derivative positive but finite and nonuniform.

## E. Falsifiability and remaining steps
- Reproduce 2021's galaxy curves using **this** truncated critical profile, its quantum normalization and *extra cutoff* parameter instead of importing 86 eV without refitting.
- Prove existence/uniqueness beyond the four discrete numerical examples, check branch turning points, energy–Casimir restrictions and radial spectral eigenvalues.
- Compare against King/fermionic King and bounded lowered Bose prior art, verifying mathematical novelty before publication.
- Include finite \(g\), coherent Bose-star exchange and actual formation kinetics — the static lowered occupation does not solve the cosmological origin problem.
- Check dependence on \(\eta\): the exact critical exponent and linear coefficient vary; finite-mass positivity depends on the chosen energy cutoff and all assumptions.

References: Guo & Rein (2007), *A Non-Variational Approach to Nonlinear Stability ... King*, https://doi.org/10.1007/s00220-007-0212-8; Kirejczyk et al. (2022), https://doi.org/10.1103/PhysRevE.105.034145; original JCAP 2021.
