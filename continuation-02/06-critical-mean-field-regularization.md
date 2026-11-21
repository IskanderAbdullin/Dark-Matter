# Manuscript 6 — A nonuniform weak-coupling limit of near-critical self-gravitating Bose halos
**New derivation and numerical checks, 2026-11-21; not peer-reviewed, not novelty-cleared.**

## Abstract
We show that neglect of repulsive Hartree–Fock interactions in the central region of a critically degenerate noncondensed Bose halo is mathematically nonuniform: the ideal model yields a finite-density radial linear cusp while **every fixed positive contact coupling** makes its leading radial density correction quadratic. The regularized region, however, scales linearly with the dimensionless interaction parameter and may be entirely unobservable for particle parameters used in Abdullin & Popov (2021). We derive the implicit equation of state, its local asymptotics, and a fiducial physical estimate, and numerically integrate the corresponding potential-based Poisson–Bose equations. This is a static phenomenological result, not a claim of dynamically realized equilibrium.

## 1. Start with the term actually present in the source
Source equations (distribBE), (nonconddens), (z), (ncEq), and (Bose-gas) have energy \(p^2/(2m)+m\Phi+2gn\), with \(g=4\pi\hbar^2 a/m\) and \(\sigma=8\pi a^2\). The original simplifies the equation by dropping \(g\) as small on halo scales. But its density equation contains the factor \( [2g/(\lambda_{\rm dB}^3kT)]\operatorname{Li}_{1/2}(z)\), **not merely** the prefactor: \(\operatorname{Li}_{1/2}(z)\to+\infty\) as \(z\to1^{-}\). The discarded term is nonuniform exactly where the ideal calculation predicts a linear central correction.

Define \(V_0^2=kT/m\), \(\nu=n_0\lambda_{\rm dB}^3\), \(R^2=V_0^2/(4\pi Gmn_0)\), \(x=n/n_0\), \(\xi=r/R\), \(y=(\Phi-\Phi_0)/V_0^2\), and \(b=2gn_0/(mV_0^2)>0\). At **central critical fugacity** \(z_0=1\), \(\nu=\nu_c=\zeta(3/2)\) and the local effective fugacity obeys
\[
z=e^{-t},\qquad t=y+b(x-1),\qquad
x=\frac{\operatorname{Li}_{3/2}(e^{-t})}{\nu_c}.
\tag{1}
\]
Equivalently \(y=t+b[1-x(t)]\) with \(t\ge0\), giving a unique nonnegative \(t\) for every \(y\ge0\); \(y(t)\) is strictly increasing. Poisson's equation becomes
\[
(\xi^2 y')'=\xi^2x(y),\qquad y(0)=y'(0)=0.
\tag{2}
\]
As a useful check, when \(b=0\), (1) reduces to the potential formulation of the original ideal halo model.

## 2. Proposition and proof
Let \(c_*=2\sqrt{\pi}/\nu_c\). For \(t\downarrow0\),
\[
1-x(t)=c_*\sqrt{t}+O(t).
\tag{3}
\]
Write \(s=1-x\). From (1), \(y=t+bs\), so \(s\sim y/b\) when \(y\to0\) with **fixed** \(b>0\): indeed \(t=O(s^2)\) and thus \(y=bs+O(s^2)\). Inverting,
\[
x(y)=1-\frac{y}{b}+O_b(y^2).
\tag{4}
\]
A finite central density yields \(y(\xi)=\xi^2/6+O(\xi^4)\) for the regularized problem. Combining (4) with Poisson gives
\[
x(\xi)=1-\frac{\xi^2}{6b}+O_b(\xi^4),\qquad
y(\xi)=\frac{\xi^2}{6}-\frac{\xi^4}{120b}+O_b(\xi^6).
\tag{5}
\]
The ideal \(b=0\) limit, in contrast, follows directly from (3):
\[
x(\xi)=1-\frac{2\sqrt\pi}{\sqrt{6}\,\nu_c}\xi+O(\xi^2)
=1-0.553979545797\,\xi+O(\xi^2)\quad(\xi\ge0).
\tag{6}
\]
Therefore \(\xi\to0\) and \(b\to0^+\) do not commute for the leading order of the central density. The finite-coupling profile is differentiable as a 3-D Cartesian scalar field at the origin; the ideal critical profile is continuous but not differentiable there.

For fixed subcritical \(z_0<1\), the analogous linear response is \(x(y)=1-[A_0/(1+bA_0)]y+O(y^2)\), \(A_0=\operatorname{Li}_{1/2}(z_0)/\nu\); the sharp critical result is the limiting \(A_0\to\infty\).

## 3. Cross-over
The quadratic and ideal square-root asymptotics switch when \(y\) is of order \((bc_*)^2\); since \(y\sim\xi^2/6\), the **order-of-magnitude** crossover is \(\xi_*\sim\sqrt6\,c_* b\). This is not a unique fitted transition radius; numerical prefactors depend on how 'crossover' is defined. At fixed \(\xi>0\), the limit \(b\to0\) recovers the ideal solution.

## 4. Fiducial physical scale, with assumptions explicit
Use the SOURCE'S published F583-4 parameters \(R=0.6\) kpc, \(V_0=44.1\) km/s, take particle mass \(m=86\) eV and illustrative upper cross section \(\sigma/m=1.25\,{\rm cm}^2/{\rm g}\) (source eqs. g, crosssec_via_g, crosssec_to_mass). Adopt identical-boson \(s\)-wave \(\sigma=8\pi a^2\). From \(n_0m=V_0^2/(4\pi GR^2)\):
- \(\rho_0=6.7648473\times10^{-24}\,{\rm g}/{\rm cm}^3\);
- \(a=8.73210\times10^{-19}\) metres;
- \(b=2.3560770\times10^{-13}\);
- rough \(\xi_*\simeq7.8313\times10^{-13}\), hence \(r_*\sim4.70\times10^{-10}\) pc (~14,500 km).
A simple classical scattering-time benchmark \(1/[\rho_0(\sigma/m)V_0]\simeq0.85\) Gyr is **not** the true Bose-enhanced kinetic equilibration time, nor proof of a common chemical potential. The cross section is a chosen illustrative bound, not a measured value for the particle. At this extreme \(b\), central regularization is mathematically exact but observationally irrelevant on usual kpc scales.

## 5. Independent computations / checks
The companion Python program constructs the critical equation of state from \(\operatorname{Li}_{3/2}(e^{-t})\), interpolates the monotone \(y(t)\), and integrates (2) in potential-flux variables. Eight numerical checks passed in a separate execution (including strict suppression of the ideal deficit by repulsion, small-radius limit, and convergence to the ideal profile away from the inner boundary layer). Representative dimensionless density at \(\xi=0.03\): \(x_{b=0}=0.9835326\), \(x_{b=0.001}=0.9844079\), \(x_{b=0.01}=0.9902481\), \(x_{b=0.05}=0.9970945\). The results are numerical **toy profiles** at artificial values of \(b\), not galaxy measurements.

## 6. Interpretative and novelty caveats
The weakness of a small-contact-coupling expansion at \(z\simeq1\) has analogues in the wider theory of interacting Bose gases. The mathematical mechanism might already exist in this literature: this note does not assert first discovery. Critical finite-temperature approximations can break down close enough to condensation; a full many-body treatment and dynamical condensate handling may modify the microscopic boundary layer. Neither a small positive \(g\) nor this static local derivation solves the cosmological free-streaming problem.

**Sources:** Abdullin & Popov, JCAP 2021 055 https://doi.org/10.1088/1475-7516/2021/11/055; Abdullin, *Dark matter halo from light bosons* (2025) https://doi.org/10.17238/issn2226-8812.2025.1.7-11; Chavanis (2025) https://doi.org/10.3389/fspas.2025.1538434; Bose-gas critical analysis https://www.phys.ens.psl.eu/~laloe/Publications/comm-bec.pdf.