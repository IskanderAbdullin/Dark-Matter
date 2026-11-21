# 14. A nonuniform infrared phase-volume exponent under arbitrarily weak Hartree–Fock repulsion

**Analytical derivation plus numerical check of the *local* coupled equation of state, 21 November 2026.** This is a subtle idealized mean-field critical endpoint, not a physical prediction about a galactic Bose-star population or a proven exact quantum many-body scaling law.

## 1. Self-consistent contact interaction modifies not only density, but the singular momentum–space geometry

Set \(\nu_c=\zeta(3/2)\), \(V_0^2=k_BT/m\), \(\xi=r/R\), \(y=[\Phi(r)-\Phi(0)]/V_0^2\), and \(b=2gn_0/(mV_0^2)>0\). In the simple repulsive Hartree–Fock description of the **noncondensed** central-critical gas (without compact Bose stars), define \(x=n/n_0\) and
\[
t(\xi)=y(\xi)+b[x(\xi)-1],\quad
x=\frac{\operatorname{Li}_{3/2}(e^{-t})}{\nu_c},\quad
f(q,\xi)=\frac{1}{\exp[q^2+t(\xi)]-1}.
\]
Earlier continuation-02 paper 06 showed that for *any fixed* \(b>0\), the central spatial density has a quadratic rather than ideal linear leading correction, although the physical boundary layer can be microscopic. The present analysis concerns the **high-occupation volume in six-dimensional phase space**.

Let \(c_*=2\sqrt\pi/\nu_c\). As \(t\downarrow0\), \(1-x=c_*\sqrt t+O(t)\). Hence \(y=t+b(1-x)\sim b c_*\sqrt t\) for any fixed \(b>0\) when \(y\to0\), so
\[
\boxed{t\sim\frac{y^2}{b^2c_*^2}.}
\]
Since regular Poisson gravity gives \(y(\xi)=\xi^2/6+o(\xi^2)\), we obtain
\[
\boxed{t(\xi)\sim\frac{\xi^4}{36b^2c_*^2}.}\tag{1}
\]
Thus the occupation singularity changes from \(f\sim[q^2+\xi^2/6]^{-1}\) at \(b=0\) to
\[
\boxed{f\sim\left[q^2+\frac{\xi^4}{36b^2c_*^2}\right]^{-1}}\qquad
(\xi\to0,\ b>0\ \mathrm{fixed}).\tag{2}
\]
The density is *spatially smoother*, yet the occupation is high across a **larger spatial region** for each fixed extremely high F. These are genuinely noncommuting limits \(b\to0\) and \(F\to\infty\).

## 2. Exact leading phase-volume power and an unexpected \(L^p\) threshold

Set \(L_F=\ln(1+1/F)\sim F^{-1}\). The superlevel condition \(f>F\) means \(q^2+\xi^4/[36b^2c_*^2]<L_F\) in the strict local asymptotic regime. Angular integration in 6D and \(s=\xi/[36b^2c_*^2 L_F]^{1/4}\) yield
\[
\boxed{
\frac{V_b(F)}{R^3(mV_0)^3}
\sim32\sqrt3\,\pi^2\,
B\!\left(\frac34,\frac52\right)
(bc_*)^{3/2}\, L_F^{9/4},
}\tag{3}
\]
where \(B(a,b)\) is Euler's beta function. By contrast, the ideal \(b=0\) critical halo has \(V_0(F)\sim4\sqrt3\pi^3 R^3(mV_0)^3 L_F^3\).

Using the layer-cake formula of paper 13, for a finite-volume cutoff compatible with this **local** center, the high-F contribution to \(\int f^p d\Gamma\) converges if and only if
\[
\boxed{p<9/4\quad (b>0\ {\rm fixed}),\qquad
p<3\quad (b=0).}\tag{4}
\]
At the equality the divergence is logarithmic. The Hartree–Fock smoothing of the macroscopic *density* does **not** mean that every high-occupation moment becomes more regular. This does not automatically imply a divergent real scattering rate or invalidate the thermodynamic EOS: the actual collision integral has momentum constraints, cancellations and coherent modes absent here.

## 3. Independent local implicit-EOS computation

Using high-precision mpmath to solve \(y=t+b[1-\operatorname{Li}_{3/2}(e^{-t})/\nu_c]\) at \(b=0.01\), the predicted limit \(t/y^2\to (bc_*)^{-2}=5430.7685\) is approached:

| \(y\) | Direct \(t/y^2\) | Direct / predicted |
|---:|---:|---:|
| 10⁻⁴ | 2815.50 | 0.5184 |
| 10⁻⁵ | 4913.08 | 0.9047 |
| 10⁻⁶ | 5372.89 | 0.98934 |
| 10⁻⁷ | 5424.91 | 0.99892 |
| 10⁻⁸ | 5430.18 | 0.999892 |

Direct beta-function angular integration independently gives high-F asymptotic slopes \(-2.99986\) for the ideal quadratic potential and \(-2.24990\) for the repulsive quartic local potential, approaching the exact exponents \(-3\) and \(-9/4\). Tests of the local asymptotics passed. This is an exact **scaling check of the leading local potential**, not a finite-radius self-consistent Hartree–Fock galactic simulation.

## 4. Why this result is almost certainly *not directly observable* at 86 eV

With the illustrative F583-4 source \(R=0.6\) kpc, \(V_0=44.1\) km/s and contact cross section \(\sigma/m=1.25\) cm²/g, the earlier result is \(b\simeq2.36\times10^{-13}\), \(c_*\simeq1.35697\). The central Hartree–Fock asymptotic requires \(y\ll(bc_*)^2\sim10^{-25}\). The corresponding occupation thresholds at which the \(F^{-9/4}\) asymptotic overtakes the ideal \(F^{-3}\) structure are **enormous**, schematically \(F\gg10^{25}\). Finite-size, HFB critical fluctuations, a coherent condensate and microscopic validity of Hartree–Fock likely matter long before such a formal continuum limit can be measured. The exponent is a mathematical sensitivity/regularization diagnostic, not a credible existing astrophysical signal or an independent particle constraint.

## 5. Next work

Repeat with a **fully specified finite-temperature HFB/GP–kinetic EOS** rather than imposing a static \(2gn\) mean field arbitrarily close to condensation, and compare to direct quantum kinetic cutoff and finite-volume regularization. Include the global tidally truncated phase-space halo and actual particle number conservation before using any \(L^p\) regularity assumptions in a nonlinear Vlasov stability theorem.

Sources: original source Abdullin & Popov JCAP 2021 055 https://doi.org/10.1088/1475-7516/2021/11/055; critical finite-temperature Bose field treatment https://doi.org/10.1103/PhysRevD.106.045025; coupled GP+Boltzmann framework https://doi.org/10.1103/PhysRevD.108.083513; our continuation-02 paper 06 and continuation-04 papers 03,10,13.
