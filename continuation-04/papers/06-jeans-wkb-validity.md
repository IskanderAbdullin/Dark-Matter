# 06 — Does the 118 pc Bose–Jeans wavelength fit inside its own high-susceptibility region?

**Date:** 2026-11-21. Original source: Abdullin & Popov (2021), eq. (Bose-gas), (inhR), (distribBE). This follows note 01 of this continuation and **qualifies** its local homogeneous Jeans analysis. The key scientific distinction is between a local static susceptibility and an eigenmode of the finite, radially stratified galaxy.

## 1. The mismatch between wavelength and local uniformity

Let \(\nu(r)=\mathrm{Li}_{3/2}(z(r))\), with \(z(r)=z_0 e^{-y(r)}\) and \(y=m[\Phi(r)-\Phi(0)]/(kT)\). The source's central dimensionless density is normalized to \(x(0)=1\); the central potential is regular, so, for a **pure noncondensed gas**,
\[
y(\xi)=\xi^2/6+O(\xi^4),\quad \xi=r/R,\quad R^2=V_0^2/(4\pi G\rho_0).
\]
Near the critical central fugacity, define \(\delta_0=-\log z_0\ll1\). The *locally inferred homogeneous* static-response factor is
\[
A(r)\equiv\frac{\operatorname{Li}_{1/2}[e^{-(\delta_0+y(r))}]}
{\operatorname{Li}_{3/2}[e^{-(\delta_0+y(r))}]}
\sim \frac{\sqrt\pi}{\zeta(3/2)\sqrt{\delta_0+y(r)}}.
\]
For \(b=0\) the **homogeneous** marginal mode wavelength at the center is
\[
\lambda_J(0)=\frac{2\pi R}{\sqrt{A(0)}}.
\]
But a local plane wave of that wavelength is only justifiable where the background varies slowly across the relevant region. Quantify the region of *large susceptibility* by a radius \(r_{1/2}\) where \(A(r_{1/2})=A(0)/2\), *not* merely by the profile's density core radius. Leading near-critical expansion gives \(y(r_{1/2})\simeq3\delta_0\), \(r_{1/2}\simeq R\sqrt{18\delta_0}\).

## 2. Numerical evaluation with the original JCAP F583-4 scale

The source's printed \(\nu=2.61\) corresponds to \(\delta_0=4.4924654222\times10^{-7}\); its \(R=600\) pc.

- Independent special-function root \(A(0)=1012.63319\).
- Exact threshold \(A(\delta_0+y)=A(0)/2\) combined with central \(y\simeq(r/R)^2/6\) gives \(r_{1/2}\simeq1.70701\) pc.
- Homogeneous \(\lambda_J(0)\simeq118.4691\) pc, **about 69.4 times the half-susceptibility radius**.
- Simpler inner transition \(\sqrt{6\delta_0}R\simeq0.9851\) pc.

So a homogeneous sinusoidal perturbation of the formally unstable central wavelength would probe a region in which the susceptibility and velocity distribution have already changed drastically. It would be **incorrect to interpret 118 pc as a demonstrated unstable wavelength of the actual spherical galaxy**. The large response arises from an extremely small near-critical central phase-space region.

| Central \(\nu\) | \(r_{1/2}\) pc, central approximation | homogeneous \(\lambda_J(0)\) pc | ratio |
|---:|---:|---:|---:|
| 2.2 | 344.75 | 1550.19 | 4.50 |
| 2.5 | 83.73 | 813.32 | 9.71 |
| 2.6 | 8.9217 | 270.36 | 30.30 |
| 2.61 | 1.7070 | 118.47 | 69.40 |
| 2.612 | 0.26957 | 47.0949 | 174.70 |

**Accuracy caveat:** for \(\nu=2.2\), \(r_{1/2}/R\simeq0.57\), and the *quadratic potential approximation* is no longer controlled. The first two rows are rough diagnostics; the high-degeneracy last three rows are quantitatively more trustworthy because their radii satisfy \(r\ll R\).

## 3. The critical-center interaction layer also fails a local WKB reading

At exact central saturation \(\delta_0=0\), the ideal susceptibility behaves \(A(r)\sim r^{-1}\). The local *homogeneous* wavelength \(\lambda_J(r)\sim R\sqrt{r/R}\), while the logarithmic gradient length of the susceptibility \(L_A=|d\log A/dr|^{-1}\sim r\). Hence
\[
\frac{\lambda_J(r)}{L_A}\sim (R/r)^{1/2}\longrightarrow\infty
\qquad (r\to0).
\]
This nonuniformity is independent of the numerical value of \(R\). A tiny positive Hartree–Fock contact term \(b\) smooths the central ideal singularity: source F583-4 illustrative \(b=2.3561\times10^{-13}\) gives a contact smoothing radius of order \(4.70\times10^{-10}\) pc, while homogeneous critical HF \(\lambda_J\simeq0.0018299\) pc—roughly \(3.89\) million times the layer width. This means the formal HF homogeneous marginal length is not automatically a resolved instability/fragmentation scale *of the stratified halo* either.

## 4. What is and isn't established

**Established mathematically:** response factor in the source BE distribution, the central scale separation and divergence of local WKB wavelength divided by the susceptibility-gradient scale. Exact numerical roots for \(\nu=2.61\).

**NOT established:** global stability or global instability of actual stationary \(f(E)\), formation of compact condensates on 118 pc or 377 AU scales, or a new astrophysical exclusion. A stable Vlasov steady state can still have a formal homogeneous Jeans instability under the 'Jeans swindle'; the local and global spectral problems are distinct.

## 5. Precise next computation

- Formulate the **linearized spherical Vlasov–Poisson integral eigenproblem** on each finite-mass isotropic lowered-Bose equilibrium (continuation-04 paper 04) with \(\ell=0,1,2\), and impose mass/energy and accessible-perturbation constraints.
- Quantify finite-cutoff sensitivity and the exactly critical unbounded \(f(E)\) endpoint, where standard differentiability assumptions may fail.
- Include separate coherent component and finite \(g\) only after the collisionless stability baseline is validated.
- Check that the same perturbation solver recovers classic King stability tests and familiar isotropic stability results (Guo & Rein 2007), not simply an uncontrolled local-growth prediction.

### Prior art and source URLs
[Abdullin & Popov JCAP 2021](https://doi.org/10.1088/1475-7516/2021/11/055); [review of self-gravitating Vlasov equilibria and stability](https://doi.org/10.1088/1361-6382/acf436); [Guo–Rein 2007 nonlinear King stability](https://doi.org/10.1007/s00220-007-0212-8); [Kirejczyk et al. 2022 noncondensed Bose density profiles](https://doi.org/10.1103/PhysRevE.105.034145). These sources support the modeling distinction; this particular numeric WKB audit is an independent calculation.
