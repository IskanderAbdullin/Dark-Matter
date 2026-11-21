# Paper 1 — Universal nonanalytic central asymptotics at critical Bose degeneracy
**Independent derivation draft, 2026-11-21.** Not peer reviewed; novelty relative to the full literature remains unchecked.
**Working title:** *A critical central-density cusp without a divergent density in noncondensed bosonic dark-matter halos*

## Abstract
We derive a local analytic diagnostic of the finite-temperature, noncondensed Bose–Einstein halo model of Abdullin & Popov (2021). At the limiting central phase-space density \(\nu=\zeta(3/2)\), the equilibrium particle density remains finite but has a universal *linear* radial correction. A subcritical halo instead has a quadratic central correction. The noncommuting limits of radius and degeneracy change the leading curvature of an idealized dark-matter rotation curve. We show the coefficient explicitly, verify it with an independent integration of the potential-based spherical Poisson equation, and specify conditions under which baryons, interactions and observational resolution erase this signature.

## 1. Statistical equilibrium and dimensionless variables
Consider a nonrelativistic ideal Bose gas in an isotropic, stationary, spherical gravitational potential. Let \(V_0^2=k_BT/m\), \(z(r)=\exp[(\mu-m\Phi(r))/(k_BT)]\), and \(\lambda_{\rm dB}=\sqrt{2\pi}\hbar/(mV_0)\). Its number density is
\[
n(r)=\lambda_{\rm dB}^{-3}\operatorname{Li}_{3/2}(z(r)).
\]
Set \(n_0=n(0)\), \(\nu=n_0\lambda_{\rm dB}^3=\operatorname{Li}_{3/2}(z_0)\), \(0<\nu\le\nu_c=\zeta(3/2)\). Define
\[
R^2=\frac{V_0^2}{4\pi Gmn_0},\quad \xi=r/R,\quad y(\xi)=\frac{m[\Phi(r)-\Phi(0)]}{k_BT}.
\]
For negligible mean-field interaction \(g\), the self-consistent Poisson–Bose problem is
\[
y''+\frac{2}{\xi}y'=\frac{\operatorname{Li}_{3/2}(z_0e^{-y})}{\nu}\equiv x(\xi),\qquad y(0)=y'(0)=0. \tag{1}
\]
This potential equation is numerically preferable to a density-only equation that contains \(\operatorname{Li}_{1/2}(z)\) in coefficients diverging as \(z\to1\).

## 2. Analytic proposition and proof
**Proposition.** If \(\nu=\nu_c=\zeta(3/2)\), then, as \(\xi\downarrow0\),
\[
x(\xi)=1-A\xi+O(\xi^2),\qquad
A=\frac{2\sqrt{\pi}}{\sqrt{6}\,\zeta(3/2)}
=0.553979545797\ldots . \tag{2}
\]
**Proof.** Spherical regularity of gravitational potential and finite central density in (1) yield \(y(\xi)=\xi^2/6+O(\xi^3)\). At \(\nu_c\), \(z_0=1\). The expansion of the polylogarithm near the Bose singular point is
\[
\operatorname{Li}_{3/2}(e^{-y})
=\zeta(3/2)+\Gamma(-1/2)y^{1/2}-\zeta(1/2)y+O(y^{3/2})
=\zeta(3/2)-2\sqrt\pi\,y^{1/2}+O(y).
\]
Insert \(y=\xi^2/6+O(\xi^3)\) and divide by \(\nu_c\), which proves (2). The central *density* remains finite; however, \(n(|\mathbf x|)\) has no Cartesian gradient at the exact critical center, because \(|\mathbf x|\) has no derivative at the origin. This is distinct from an NFW \(n\propto r^{-1}\) divergent-density cusp. \(\square\)

**Subcritical result.** For \(0<\nu<\nu_c\), \(z_0<1\) and
\[
x(\xi)=1-\frac{\operatorname{Li}_{1/2}(z_0)}{6\nu}\xi^2+O(\xi^4). \tag{3}
\]
The coefficient diverges as \(\nu\to\nu_c^{-}\): one cannot exchange the \(\xi\to0\) and \(\nu\to\nu_c\) limits. With \(\delta=-\log z_0\ll1\), the crossover between quadratic and effectively linear behavior occurs when \(y\sim\delta\), or
\[
\xi_{\rm cross}\sim \sqrt{6\delta}\sim
\frac{\sqrt6}{2\sqrt\pi}(\nu_c-\nu). \tag{4}
\]
Equation (4) is leading asymptotics, not a calibrated observable radius.

## 3. Effect on circular velocity
Gauss's law gives
\[
u^2(\xi)\equiv\frac{V^2(r)}{V_0^2}=\frac1\xi\int_0^\xi x(s)s^2\,ds=\xi y'(\xi). \tag{5}
\]
Using (2),
\[
u^2(\xi)=\frac{\xi^2}{3}-\frac{A\xi^3}{4}+O(\xi^4),
\quad
u(\xi)=\frac{\xi}{\sqrt3}\left[1-\frac{3A}{8}\xi+O(\xi^2)\right]. \tag{6}
\]
The subcritical case (3) has its first deviation from solid-body rise at relative order \(\xi^2\), rather than \(\xi\). Thus, *in the ideal dark-matter-only limit*, the first derivative of \(V(r)/r\) with respect to \(r\) is a criticality diagnostic. It is **not** a claim that present galaxy observations have the radial resolution to measure it.

## 4. Independent numerical check
Integrate (1) as \(y'=q/\xi^2,\ q'=\xi^2x\), starting at a small \(\epsilon\) with \(y(\epsilon)=\epsilon^2/6\), \(q(\epsilon)=\epsilon^3/3\). In our standalone calculation the analytic critical coefficient was \(0.5539795458\), and the finite-radius numerical estimate at \(\xi=10^{-5}\) was \(0.5539786370\). Independent profiles were computed for \(\nu=0.03, 0.5, 1.5, 2.2,\nu_c\). The numerical test verifies the internal mathematics of (1); it is not a cosmological or observational test.

## 5. Physical and observational discriminants
A repulsive mean-field potential \(2gn\) changes the fugacity argument and may regularize the divergence. A finite condensate fraction would also change the inner potential. It is necessary to solve the coupled self-consistent system, not extrapolate (2) into the condensed phase. Observationally, beam smearing, pressure support, non-circular motion, inclination and luminous matter can mimic or hide a weak \(r^2\) term in \(V(r)\). A meaningful analysis must simultaneously fit baryons and data covariances, and compare critical and subcritical profiles with the same number of effective degrees of freedom.

## 6. Falsifiable predictions and next calculations
1. The potential-based ideal-gas solution must display the analytically specified \(A\) at exact criticality, with residual \(O(\xi^2)\).
2. At every fixed \(\nu<\nu_c\), the density derivative at the origin vanishes; the onset of an approximately linear radial region scales with (4) as \(\nu\uparrow\nu_c\).
3. A positive mean-field \(g\), a composite Bose-star core or finite-resolution convolution may destroy distinguishability; map this parameter domain quantitatively.
4. If baryon-corrected sub-kpc rotation curves consistently disallow the predicted critical curvature *and* favor \(\nu\approx\nu_c\) from global fits, the ideal two-component model is internally disfavored.

## Limitations and references
This is an analytical result for a **stationary isotropic equilibrium Bose gas**, not an origin theory for the equilibrium distribution. Mathematical originality must be checked against earlier Bose-potential/critical-core analyses before publication. Starting points: [Abdullin & Popov (2021)](https://doi.org/10.1088/1475-7516/2021/11/055); [Chavanis (2025), review](https://doi.org/10.3389/fspas.2025.1538434); [Chavanis (2019), predictive core–halo model](https://arxiv.org/search/?query=Chavanis+2019+predictive+model+BEC+dark+matter+halos&searchtype=all).


## Appendix A — critical asymptotics in the presence of a compact-object background (new independent result)
Now include a *smooth, isothermal spatial density* of compact Bose stars as in the two-component phenomenology, \(x_{\rm BS}=\alpha e^{-\gamma y}\). The gravitational potential obeys \(y''+2y'/\xi=x_n+\alpha e^{-\gamma y}\) and therefore
\[
y(\xi)=\frac{1+\alpha}{6}\xi^2+O(\xi^3)\quad(\xi\downarrow0).
\]
At \(\nu=\zeta(3/2)\) the **smooth noncondensed** density becomes
\[
x_n(\xi)=1-A\sqrt{1+\alpha}\,\xi+O(\xi^2),
\quad A=\frac{2\sqrt\pi}{\sqrt6\zeta(3/2)}.
\]
The Bose-star component, unlike its noncondensed companion, begins quadratically:
\[
x_{\rm BS}(\xi)=\alpha-\frac{\alpha\gamma(1+\alpha)}6\xi^2+O(\xi^3).
\]
Consequently, total density normalized to its actual central value obeys
\[
\frac{x_n+x_{\rm BS}}{1+\alpha}
=1-\frac{A}{\sqrt{1+\alpha}}\,\xi+O(\xi^2).
\]
**Interpretation:** stronger central condensate gravity sharpens the *unnormalized* noncondensed cusp coefficient (\(\propto\sqrt{1+\alpha}\)), but dilutes the *total normalized* relative slope (\(\propto1/\sqrt{1+\alpha}\)). The concentration parameter \(\gamma\) enters at the next order, not the leading slope, if the Bose-star density is analytic at the origin. This follows specifically from the phenomenological isothermal Bose-star spatial distribution; it would not remain true for a singular stellar density or a point-like central object.

Independent integration for \(\alpha=\gamma=12.8\), \(\nu=\zeta(3/2)\), at \(\xi=10^{-5}\):
predicted one-sided noncondensed slope \(A\sqrt{13.8}=2.0579426748\); directly measured \(2.0579299027\), about \(6.2\times10^{-6}\) relative difference. The agreement is a local numerical unit test, not an observational detection.

The two-component demonstration in \`code/two_component_demo.py\` also shows why a dense condensate population can increase inner speeds without increasing them at every radius: the fixed central noncondensed density undergoes a self-consistent redistribution in the deeper total potential. Compare the nondimensional speeds before evaluating any specific galaxy.


## Appendix B — nonzero repulsive mean-field coupling removes the strict ideal critical cusp
A new follow-up [derivation, calculations and limitations](../continuation-02/06-critical-mean-field-regularization.md) retains the \(2gn\) term that the source neglects for the noncondensed halo. Define \(b=2gn_0/(mV_0^2)>0\). At exactly critical central fugacity the interacting equation of state is \(x=\operatorname{Li}_{3/2}(e^{-t})/\zeta(3/2)\) with \(y=t+b(1-x)\) and Poisson equation \((\xi^2 y')'=\xi^2x\). For fixed \(b>0\), \(x=1-y/b+O_b(y^2)=1-\xi^2/(6b)+O_b(\xi^4)\). This is a smooth central profile, whereas the ideal model's \(b=0\) has the one-sided linear term derived in Section 2. The limits \(b\to0\) and \(\xi\to0\) are nonuniform; the transition region is \(\xi=O(b)\), and can be microscopic for the cross-section range used in the original paper. For the illustrative F583-4 parameters, \(b\simeq2.356\times10^{-13}\) and the nominal crossover is \(r\sim4.7\times10^{-10}\) pc. This does not invalidate the ideal profile on *observable* kiloparsec radii. The new calculation is not a proof of kinetic equilibration, dynamic stability, or originality against all Bose-gas literature.
