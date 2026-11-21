# 01. Near-critical Bose gas: a collisionless *local* Jeans instability diagnostic

**Status:** independent theoretical calculation plus numerical integrations. Not peer-reviewed; global halo stability *not* established. Original model source: Abdullin & Popov (2021), eqs. (distribBE), (dispV0), (nonconddens), (Bose-gas), (inhR).

## A. Derivation

Let \(V_0^2=kT/m\), central density \(\rho_0=mn_0\), \(R^2=V_0^2/(4\pi G\rho_0)\) and \(\nu=\mathrm{Li}_{3/2}(z_0)\in(0,\zeta(3/2)]\). The positive dimensionless BE occupation is \(f_0(v)=[z_0^{-1}\exp(v^2/(2V_0^2))-1]^{-1}\). The usual uniform Jeans-background subtraction is applied; this is *not* a true equilibrium of infinite homogeneous Newtonian gravity without that prescription.

Linearize Vlasov–Poisson about a homogeneous isotropic distribution, first ignoring collisions and any coherent component:
\[
\partial_t F+\mathbf v\cdot\nabla_xF-\nabla\Phi\cdot\nabla_vF=0,\quad
\nabla^2\delta\Phi=4\pi G\delta\rho.
\]
At the **marginal zero frequency** \(\omega=0\), the density response to a perturbation \(\delta U\) of specific potential energy is
\[
\delta\rho=-\chi_0\delta U,\quad
\chi_0=-\int d^3v \frac{\partial_{v_\parallel}F_0}{v_\parallel}
=\frac{\rho_0}{V_0^2}\,
\frac{\mathrm{Li}_{1/2}(z_0)}{\mathrm{Li}_{3/2}(z_0)}.
\]
Here \(F_0\) is normalized to **mass density per velocity volume**, while \(f_0\) is the occupancy; distinguishing them prevents a false \(m\)-factor. This identity follows directly from \(\partial_{v_\parallel}f_0=-(v_\parallel/V_0^2)f_0(1+f_0)\).

Set \(A=\mathrm{Li}_{1/2}(z_0)/\nu\). For a repulsive contact Hartree–Fock shift \(2g\,n\), the specific mean-field potential is \((2g/m^2)\rho\), hence \(\delta U=\delta\Phi+(2g/m^2)\delta\rho\). Define \(b=2g n_0/(mV_0^2)\). Eliminating \(\delta\Phi\) and \(\delta\rho\) gives
\[
\boxed{(k_JR)^2=\frac{A}{1+bA}},\qquad
\boxed{\lambda_J=2\pi R\sqrt{\frac1A+b}}. \tag{1}
\]
This matches the **isothermal zero-frequency susceptibility** \(\lambda_{J,T}=2\pi c_T/\sqrt{4\pi G\rho_0}\) with
\[
c_T^2=V_0^2\left[\frac{\mathrm{Li}_{3/2}(z_0)}
{\mathrm{Li}_{1/2}(z_0)}+b\right].
\]
This equality of static response scales does **not** identify the full time-dependent collisionless dispersion relation with a fluid sound-wave law.

Classical \(z_0\ll1\): \(A\to1\), \(\lambda_J\to2\pi R\). Approaching ideal criticality \(z_0\uparrow1,b=0\): \(A\to\infty\), \(\lambda_J\to0\) in the mathematical homogeneous approximation, reflecting an infrared-divergent static susceptibility. At fixed \(b>0\), \(\lambda_J\to2\pi R\sqrt b\), *not zero*. Interactions, finite volume, quantum coherence, distribution regularity and global halo stratification are necessary before inferring any physical small-scale collapse.

## B. Results at the original galaxy F583-4 scales

Illustrative normalization \((R,V_0)=(0.6\,{\rm kpc},44.1\,{\rm km/s})\); the printed \(\nu=2.61\) is not exactly saturated. Independent numerical evaluation:

| Central \(\nu\) | \(A(z_0)\) | Ideal local \(\lambda_J\) (pc) |
|---:|---:|---:|
| 0.03 | 1.01073 | 3749.85 |
| 0.50 | 1.21848 | 3415.24 |
| 1.00 | 1.57212 | 3006.69 |
| 1.50 | 2.24317 | 2517.10 |
| 2.20 | 5.91419 | 1550.19 |
| **2.61** | **1012.633** | **118.469** |
| \(\zeta(3/2)\), \(g=0\) | divergent | zero mathematical limit |

For **exact** critical fugacity and illustrative \(\sigma/m=1.25\) cm²/g, the Hartree–Fock parameter is \(b=2.3561\times10^{-13}\), so \(\lambda_J=0.0018299\) pc (~377 AU). This is a limit of a **homogeneous** response at fixed density and coupling; it is *not* a length at which quantum stars must form, or a stability verdict for the source halo. For \(\nu=2.61\), \(bA\approx2.4\times10^{-10}\): the same illustrative contact term barely changes the 118-pc marginal wavelength.

The independent direct \(q^2 f(1+f)\) momentum quadrature at \(\nu=2.61\) gave \(A=1012.6331898322308\), versus \(1012.633189832231\) via \(\operatorname{Li}_{1/2}/\operatorname{Li}_{3/2}\), within roundoff.

## C. Full imaginary-axis kinetic dispersion for a homogeneous comparison

Seek an exponentially growing mode \(\omega=i\gamma,\gamma>0\), with \(s=\gamma/(\sqrt2\,k V_0)\). The angularly integrated response is
\[
{\cal A}(s)=
\frac{\displaystyle\int_0^\infty dq\,q^2 f(q)[1+f(q)]
\left[1-\frac{s}{q}\arctan\frac{q}{s}\right]}
{\displaystyle\int_0^\infty dq\,q^2 f(q)},\quad
f(q)=[e^{\delta+q^2}-1]^{-1},\ \delta=-\log z_0.
\tag{2}
\]
The collisionless growing-mode dispersion condition is \((kR)^2={\cal A}(s)\) when \(g=0\), or \((kR)^2={\cal A}(s)/[1+b{\cal A}(s)]\) for the same instantaneous repulsive mean field. Eq. (2) obeys \({\cal A}(0)=A\) and \({\cal A}(s\to\infty)\to1/(2s^2)\), so at \(kR\to0\), \(\gamma R/V_0\to1\), a further independent check.

For \(\nu=2.61,b=0\), numerical root-finding and quadrature yielded:

| \(kR\) | \(\lambda/R\) | \(\gamma R/V_0\) | \(e\)-fold time at printed F583-4 scales |
|---:|---:|---:|---:|
| 5 | 1.257 | 0.18320 | 72.6 Myr |
| 10 | 0.628 | 0.08603 | 154.6 Myr |
| 15 | 0.419 | 0.04965 | 267.9 Myr |
| 20 | 0.314 | 0.02899 | 458.8 Myr |
| 25 | 0.251 | 0.01468 | 906 Myr |
| 30 | 0.209 | 0.003557 | 3.74 Gyr |
| \(>\sqrt A=31.8219\) | <0.1974 | no \(\gamma>0\) found | — |

**Warning about applicability:** This is a strictly local **uniform-medium** dispersion, not the eigenvalue spectrum of a finite spherical galaxy. Modes with \(\lambda\gtrsim R\) are *especially* outside local WKB justification; even \(\lambda\sim0.2R\) has non-negligible radial stratification. Unstable homogeneous modes do not prove unstable finite \(f(E)\)-halos; boundary conditions, conserved integrals, anisotropy, baryons, a Bose-star component and wave/quantum kinetics can alter conclusions. Growth times are diagnostics, not predictions of dwarf counts or condensate nucleation.

**Testable next step:** derive and numerically solve the nonlocal \(l=0,1,2\) linearized spherical Vlasov–Poisson eigenproblem about the *actual* source \(F(E)\), and repeat for finite-mass lowered Bose and interacting two-component versions. Compare globally conserved mass, energy and Casimir perturbation classes. Do not equate marginal homogeneous \(\lambda_J\) to a galactic core size.

**Prior art / citations:** standard kinetic Jeans framework and warnings about homogeneous backgrounds: e.g. *Stability and instability results for equilibria of a self-gravitating collisionless gas — a review* (2023), https://doi.org/10.1088/1361-6382/acf436 ; Kirejczyk et al. (2022) https://doi.org/10.1103/PhysRevE.105.034145; Chavanis (2025) https://doi.org/10.3389/fspas.2025.1538434 ; Khlebnikov (2000) *Short-scale gravitational instability in a disordered Bose gas*, https://doi.org/10.1103/PhysRevD.62.043519 (different nonequilibrium/wave model). Mathematical identity here is an **application** of known kinetic principles to source JCAP; novelty across all Bose kinetic studies unverified.
