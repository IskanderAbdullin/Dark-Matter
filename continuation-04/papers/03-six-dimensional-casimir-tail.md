# 03. Six-dimensional high-occupancy Casimirs: an unexpected match in the leading power and a mismatch at the next order

**Independent analytic and numerical theoretical note, 2026-11-21.** It does not prove physical relaxation or a collisionless origin for Abdullin & Popov's 2021 equilibrium; novelty against kinetic literature needs further assessment.

## A. Why comparing only the momentum-space singularity is misleading

A relativistically decoupled, homogeneous, zero-chemical-potential Bose relic has \(f_{\rm ini}(q)=(e^{q/p_0}-1)^{-1}\sim p_0/q\) at \(q\to0\), whereas a *local nonrelativistic critical* Bose distribution has \(f_{\rm fin}(p,r=0)\sim 2m^2V_0^2/p^2\). Naively their large-occupation momentum-space phase-volume tails scale differently: \(F^{-3}\) versus \(F^{-3/2}\). **But a galaxy's critical fugacity occurs at a spatial point, not throughout a homogeneous three-dimensional volume.** The correct invariant of collisionless phase flow is a **six-dimensional** volume in *both positions and momenta*.

For a dimensionless occupation \(f\) define the superlevel phase-space volume
\[
{\cal V}_f(F)=\int d^3x\,d^3p\;\mathbf{1}_{f(\mathbf{x},\mathbf{p})>F}, \quad F>0.
\tag{1}
\]
The optional common factor \((2\pi\hbar)^{-3}\) counts quantum cells. Under an **exact, smooth, closed Hamiltonian/Vlasov map** \(f(t)\) is constant along trajectories and the canonical volume is preserved, so (1) must be constant for every \(F\). This is a strong, fine-grained Casimir statement. Observed or coarse-grained haloes need not conserve (1) for the **coarse-grained** \(f\).

## B. Initial relativistic freeze-out spectrum

For one homogeneous Lagrangian comoving patch \(V_L\), \(f_i(q)=[e^{q/p_0}-1]^{-1}\) at \(\mu=0\); with present canonical comoving momentum \(q\), \(\log(1+1/F)\equiv L_F\). Then \(f_i>F\iff 0<q<p_0L_F\), exactly:
\[
\boxed{{\cal V}_i(F)=V_L\,\frac{4\pi}{3}\,p_0^3\,L_F^3.}\tag{2}
\]

## C. Final source critical spherical BE halo (six dimensions)

At \(z_0=1\), \(f_f(r,p)=[e^{u^2+y(\xi)}-1]^{-1}\), \(u=p/(\sqrt2mV_0)\), \(\xi=r/R\). Poisson with finite central density implies \(y(\xi)=\xi^2/6-A_c\xi^3/12+O(\xi^4)\), where \(A_c=2\sqrt\pi/(\sqrt6\,\zeta(3/2))\). This differs from an exactly homogeneous *critical* medium. The \(f_f>F\) six-dimensional region obeys
\[
u^2+y(\xi)<L_F.
\]
Its phase-space volume is, by direct angular integration,
\[
{\cal V}_f(F)=\frac{16\pi^2}{3}(\sqrt2mV_0)^3R^3
\int_0^{\xi_F}\xi^2[L_F-y(\xi)]^{3/2}d\xi,\quad
y(\xi_F)=L_F. \tag{3}
\]
At \(F\to\infty\), using \(y\simeq\xi^2/6\), this region is a 6-dimensional ellipsoid, and
\[
\boxed{{\cal V}_f(F)\sim 4\sqrt3\,\pi^3 R^3(mV_0)^3 L_F^3
\propto F^{-3}.}\tag{4}
\]
**The leading high-occupation exponent matches the primordial relativistic relic's \(-3\)**. A claim that the spectra cannot be collisionlessly related *solely because* the momentum-space singularities are \(1/q\) and \(1/p^2\) would be wrong. Their full six-dimensional geometry matters.

## D. A subleading term proves the two ideal spectra are not *exactly* level-volume identical in a closed patch

Write \(\xi=\sqrt{6L_F}\,t\). The cubic correction to \(y\) gives \(y/L_F=t^2-c\sqrt{L_F}t^3+\cdots\), \(c=A_c\sqrt6/2=\sqrt\pi/\zeta(3/2)\). Using \(\int_0^1t^2(1-t^2)^{3/2}dt=\pi/32\), \(\int_0^1t^5(1-t^2)^{1/2}dt=8/105\), expansion of (3) yields
\[
\boxed{
{\cal V}_f(F)=4\sqrt3\pi^3 R^3(mV_0)^3 L_F^3
\left[1+B\sqrt{L_F}+O(L_F)\right],
\quad B=\frac{128}{35\sqrt\pi\,\zeta(3/2)}=0.789826\ldots
} \tag{5}
\]
In (2) the ratio \({\cal V}_i(F)/L_F^3\) is **exactly constant** for a homogeneous zero-\(\mu\) relativistic Bose relic. In (5), the normalized ratio varies with \(F\) at order \(F^{-1/2}\). Therefore no single choice of fixed \(V_L\) makes \({\cal V}_i(F)={\cal V}_f(F)\) at all high \(F\), **if**: initial distribution is exactly homogeneous thermal-relativistic BE, final one is exactly the fine-grained ideal source critical spherical BE profile, the same closed collection of particles is followed, and no collisions, selective escape or coarse-graining occur. This is a **conditional exact-fine-grained mismatch**, not a rejection of any coarse-grained galactic BE fit. Allowing primordial nonthermal \(f(q)\), inhomogeneity, particle exchange, collisions, coherent condensate or incomplete mixing can evade the premise.

## E. Independent numerical verification

Poisson was independently integrated for the exact-critical density \(x(y)=\mathrm{Li}_{3/2}(e^{-y})/\zeta(3/2)\) using a controlled small-y polylog expansion and a root event at \(y=L_F\). The **full six-dimensional** integral (3) was evaluated numerically.

| Occupancy threshold F | Measured \({\cal V}_f/(C L_F^3)\) | First-order formula \(1+B\sqrt{L_F}\) |
|---:|---:|---:|
| 10 | 1.296867 | 1.243838 |
| 100 | 1.083702 | 1.078786 |
| 1000 | 1.025446 | 1.024970 |
| 10,000 | 1.007945 | 1.007898 |
| 100,000 | 1.002503 | 1.002498 |

Here \(C=4\sqrt3\pi^3R^3(mV_0)^3\). At F=100000 the relative first-order error is ~\(5\times10^{-6}\). The numerical slope on high-F points is \(-3.00437\). The numerical experiment supports (4)–(5) but neither measures phase-space mixing in a real halo nor verifies a future collision solver.

## F. A formal Lagrangian-mass anchor (interesting but extremely assumption-dependent)

If one matches *only* the leading asymptotic coefficients in (2) and (4), with the same scalar species:
\[
V_L=3\sqrt3\,\pi^2 R^3\left(\frac{mV_0}{p_0}\right)^3.
\]
For a **pure, one-real-scalar zero-\(\mu\)** relic that constitutes the relevant primordial DM population, \(n_{b0}=\zeta(3)p_0^3/(\pi^2\hbar^3)\), so the formal patch mass is
\[
\begin{aligned}
M_L&=m n_{b0}V_L
 =3\sqrt3\,\zeta(3)\frac{m^4V_0^3R^3}{\hbar^3}\\
 &=\frac{3\sqrt3\,\zeta(3)\sqrt{\pi/2}}{\nu_c}
 \frac{V_0^2R}{G}
 \simeq2.996617\,\frac{V_0^2R}{G},\quad(\nu=\nu_c).
\end{aligned}\tag{6}
\]
The formal mass is nearly three times the halo's natural mass scale \(M_0=V_0^2R/G\), and the arbitrary relic momentum temperature cancels. At \(R=0.6\) kpc and \(V_0=44.1\) km/s: \(M_0\simeq2.7130\times10^8 M_\odot\), and \(M_L\simeq8.1299\times10^8 M_\odot\) for exactly critical \(m=86.20496\) eV, with an illustrative thermal \(p_0\simeq9.75\times10^{-5}\) eV/c corresponding to a comoving Lagrangian radius ~179.8 kpc.

**Do not use (6) to infer a cosmological halo mass.** It matches a leading tail coefficient under highly restrictive closed-patch assumptions; the next-order Casimir mismatch already shows exact pure collisionless mapping of the whole distributions is not achieved. Infinite-mass original halo, mixing with neighbouring Lagrangian patches, baryons and cosmic nonuniformity make the patch matching illustrative, not an observed mass prediction.

## G. Next scientific step

Derive the phase-space *level-volume function* \({\cal V}_i(F)\) of a physically produced cold Bose \(f(q)\) in the Yin/Sakurai–Yin mechanism, including a finite primordial cutoff and any coherent component. Compare to self-consistent finite-mass Bose halo \({\cal V}_f(F)\) across the entire F range. If incompatible under fine-grained collisionless evolution, calculate the **minimum required collision/coherent exchange** or permitted mixing; do not simply assume the late BE shape emerges from violent relaxation.

References: JCAP source (2021) https://doi.org/10.1088/1475-7516/2021/11/055; collisionless invariants and relaxation Chavanis (2022) https://doi.org/10.1016/j.physa.2022.128089; 2023 hybrid quantum kinetic evolution https://doi.org/10.1103/PhysRevD.108.083513; 2025 stochastic formulation https://doi.org/10.1103/PhysRevD.111.023505. Newness relative to all prior kinetic/Casimir asymptotics **not established**.
