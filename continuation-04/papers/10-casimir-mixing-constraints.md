# 10. Phase-space mixing constraints on the *formation* of finite-mass critical Bose halos

**Independent theory and numerical experiment (21 Nov 2026).** Based on the source Abdullin–Popov 2021 stationary Bose model, our *positive finite-mass lowered Bose* extension (papers 04,07,09), and the mixing theorems of Tremaine–Hénon–Lynden-Bell (1986) and Dehnen (2005). This is a **conditional necessary-condition analysis**, NOT a general exclusion of an 86-eV bosonic particle or proof that the distribution is dynamically realized.

## 1. Fine-grained Casimirs vs coarse-grained mixing

For an isotropic occupation \(f\ge0\) define the six-dimensional superlevel volume \(V_f(F)=\int d^3x\,d^3p\,\mathbf1_{f>F}\). A closed, exactly collisionless Hamiltonian evolution conserves this at **every** \(F\). A coarse-grained final \(\bar f\) need not conserve it; a more appropriate necessary condition (Dehnen's excess-mass mixing theorem, in its domain of validity) uses
\[
D_f(F)=\int [f(\mathbf x,\mathbf p)-F]_+\,d^3x\,d^3p
=\int_F^\infty V_f(u)\,du .
\]
Under phase-space mixing of a closed initial population,
\[
\boxed{D_{\rm final}(F)\le D_{\rm initial}(F)\quad\text{for all }F>0.}\tag{1}
\]
This is a necessary comparison of **the entire distributions**, not just their central occupation or a chosen six-dimensional volume. The common factor \((2\pi\hbar)^{-3}\) cancels for a single particle species. It requires coherent choices of canonical coordinates; in homogeneous cosmology the comoving pair \((\mathbf x_{\rm com},\mathbf q=a\mathbf p_{\rm phys})\) has the same phase-volume element as the physical pair.

## 2. Explicit initial and final models

**Initial:** one species of collisionlessly decoupled, zero-chemical-potential *relativistic* Bose particles, exactly homogeneous in a finite comoving Lagrangian patch \(V_L\):
\[
f_i(q)=(e^{q/p_0}-1)^{-1}.
\]
Then, with \(L_F=\ln(1+1/F)\),
\[
V_i(F)=\frac{4\pi}{3}V_Lp_0^3L_F^3.
\]
The patch's initial particle number is \(N_i=V_L\,\zeta(3)p_0^3/(\pi^2\hbar^3)\).

**Final:** our exactly central-critical (\(\psi_0=\eta\)) lowered Bose halo, bounded in both physical radius \(r_t=R\xi_t\) and particle energy, with \(f_f=[e^{q_p^2+y(\xi)}-1]^{-1}-c_\eta\) inside the escape boundary, \(q_p=p/(\sqrt2mV_0)\), \(c_\eta=(e^\eta-1)^{-1}\), \(y(\xi)=\eta-\psi(\xi)\), and \(f_f=0\) outside. All particles of the closed initial patch are assumed retained and no coherent condensate is added for this test. The actual *number* normalization is \(\nu_\eta=4I_c(\eta)/\sqrt\pi\), and total final mass \(M_f=Q_t\,V_0^2R/G\). Both quantities were calculated directly from the self-consistent Poisson–Bose solution (paper 09).

These are intentionally **strong, falsifiable hypotheses**, not claimed to follow from cosmological galaxy formation.

## 3. Exact leading high-occupation necessary bound

For \(F\to\infty\), the six-dimensional final volume has the same leading \(L_F^3\) power as the initial distribution because the critical point is localized in *both* space and momentum. With \(C_0=4\sqrt3\,\pi^3\),
\[
V_f(F)=C_0R^3(mV_0)^3 L_{F+c_\eta}^3
[1+B_\eta\sqrt{L_{F+c_\eta}}+O(L_F)],
\quad
\boxed{B_\eta=\frac{32}{35 I_c(\eta)}>0.} \tag{2}
\]
The \(F^{-7/2}\) correction follows from the critical linear density slope and a six-dimensional ellipsoidal integral. At \(\eta=4\), \(I_c=1.087257505\), hence \(B_\eta=0.8409100\). A **new high-accuracy independent** solution of the central Poisson ODE, with momentum-quadrature and adaptive integration, gives at \(F=10^5\) normalized volume ratio **1.002664486** versus first-order **1.002659184**, absolute discrepancy \(5.30\times10^{-6}\); at \(F=10^6\) discrepancy \(5.28\times10^{-7}\). An initial lower-precision ODE produced spurious \(10^{-4}\) residuals: these were traced to numerical accuracy and replaced, not reported as physical corrections. The exact leading exponent remains \(-3\).

If the initial mass is **exactly equal** to the finite halo mass (no escape), equality \(N_i=M_f/m\) fixes \(V_L\). The ratio of leading initial to final phase-volume coefficients is
\[
\boxed{s_\eta=\frac{Q_t(\eta)\nu_\eta}{3\sqrt3\,\zeta(3)\sqrt{\pi/2}}}. \tag{3}
\]
For the final mixed state not to exceed the initial **high-\(F\)** excess mass, one needs \(s_\eta\ge1\); if exactly equal, the positive \(B_\eta F^{-1/2}\) term still creates an excess, so strict inequality is needed in that idealized limit. The constant in the denominator is **7.8282889486**. In the seven-point family:
| cutoff \(\eta\) | \(Q_t\nu_\eta/7.8283\) | High-F test |
|---:|---:|---|
| 0.5 | 0.1388 | Fails |
| 1 | 0.5513 | Fails |
| 2 | 2.1533 | Passes *this necessary check* |
| 3 | 4.7166 | Passes *this necessary check* |
| 4 | 8.2613 | Passes *this necessary check* |
| 6 | 20.0006 | Passes *this necessary check* |

Independent numerical threshold within this specific family: \(s_\eta=1\) at \(\eta\simeq1.35167\) (root of the computed equilibrium family). This is **not** a physical lower bound on the real galaxy's tidal cutoff, because assumptions of a thermal homogeneous progenitor and all-mass retention are unproven. For an **open** initial patch with mass \(M_i\ge M_f\), the high-F necessary mass factor is \(M_i/M_f\ge1/s_\eta\) when \(s_\eta<1\). For \(\eta=1\), this means **at least 44.9%** of a sufficiently large progenitor's mass must be absent from the final equilibrium under this high-F-only consistency test, if all other assumptions hold.

## 4. Stronger finite-F excess-mass tests

The criterion (1) is stronger than leading high-F volume inequality. For fixed \(\eta\), direct quadrature evaluates
\[
\frac{D_i(F)}{C_0R^3(mV_0)^3}
=s_\eta\int_0^{L_F}
\frac{L^3e^L}{(e^L-1)^2}\,dL,
\]
\[
\frac{D_f(F)}{C_0R^3(mV_0)^3}
=\int_0^{L_{F+c_\eta}}\frac{L^3 e^L}{(e^L-1)^2}
\frac{V_f(L)}{C_0 R^3(mV_0)^3 L^3}\,dL,
\]
where \(V_f(L)\) is obtained by integrating the *actual numerical Poisson potential* over the \(q_p^2+y<L\) ellipsoid, not by inserting the leading asymptotic alone.

Examples with no-escape mass equality:
| \(\eta\) | \(D_f/D_i\) at \(F=0.001\) | at \(F=1\) | at \(F=100\) |
|---:|---:|---:|---:|
| 1 | 1.0684 | **2.7101** | 2.0115 |
| 2 | 1.0637 | 0.8667 | 0.5037 |
| 4 | **1.0158** | 0.2306 | 0.1296 |
| 6 | 0.8147 | 0.0938 | 0.0534 |

Thus even the \(\eta=4\) equilibrium passing the high-occupation **leading** bound violates the full coarse-grained mixing condition at \(F=10^{-3}\) under **closed thermal-patch** assumptions. This does not establish a minimum physically admissible \(\eta\), for an additional reason below.

With a progenitor patch allowed to lose particles, the sampled maximum \(D_f/D_i\) is a **lower bound** on the requisite factor by which the primordial patch particle number must exceed that of the final halo. For example at \(\eta=4\), \(D_f/D_i=1.0158\) at \(F=0.001\) implies at least \(\sim1.56\%\) mass loss under the idealized one-population setup, not a predicted real evaporation fraction. At \(\eta=1\), \(D_f/D_i=2.7101\) at \(F=1\) instead implies at least \(\sim63.1\%\) mass not retained. These are conditional, **not observed** escaping masses.

## 5. A general support theorem: exact sharp cutoff cannot arise by closed mixing from an unbounded thermal momentum tail

For *every finite \(\eta\)* the final phase-space support has finite six-dimensional volume:
\[
V_f(0^+)=\frac{16\pi^2}{3}(\sqrt2mV_0)^3R^3
\int_0^{\xi_t}\xi^2\psi(\xi)^{3/2}d\xi<\infty.
\]
For an initially relativistic Bose distribution, \(f_i(q)>0\) for arbitrarily large \(q\), so \(V_i(F)\sim C_i[\ln(1/F)]^3\to\infty\) as \(F\downarrow0\) in any nonzero finite spatial patch. With **exactly equal initial and final total particle number**, \(D(0)=N\), hence
\[
D_f(F)=N-FV_f(0^+)+o(F),\qquad
D_i(F)=N-C_iF\{[\ln(1/F)]^3+O([\ln(1/F)]^2)\}.
\]
Consequently
\[
\boxed{D_f(F)>D_i(F)\ \text{for sufficiently small }F>0}
\]
for any finite \(\eta\). This conditional no-go is fundamentally a **support/truncation mismatch**, not a special feature of 86-eV Bose microphysics. A realistic tidal halo has an escaping population and never has a mathematically exact step-function cutoff, so the theorem tells us to model a *nonzero escaping phase-space tail* or an open system, **not to reject the particle theory**.

Numerical consistency: at \(\eta=4,F=10^{-8}\), \(D_f/D_i\approx1.00000908\), tending to 1 from above as \(F\to0\). At \(\eta=6\) the violation shifts to far smaller \(F\), and a finite sampled threshold grid can miss it. Accordingly no universal physical cutoff bound is declared.

## 6. Precise research consequence

The finite-mass lowered BE profile is an excellent stationary idealization and passes hydrostatic/virial checks, yet **cannot be the *exact* all-particle-retained fine- or coarse-grained collisionless image** of a homogeneous untruncated zero-\(\mu\) thermal relic patch. A physically complete theory must include (i) an actual finite primordial spectrum, (ii) escaping or accreted particles and an explicitly open halo, (iii) collisional Bose \(C_{22},C_{12}\) and/or coherent condensate exchange, or (iv) a soft rather than hard tidal tail. The Dehnen inequality then has to be recalculated for that physically specified distribution and mass balance.

**Novelty and scope:** This is our original *application and numerical comparison* to this family of Bose halos, not a new mixing theorem. The specific analytic and numeric cutoffs have not been checked against every prior theoretical model. The directly applicable prior results are:
- Tremaine, Hénon & Lynden-Bell, *H-functions and mixing in violent relaxation*, MNRAS 219 (1986) 285, https://doi.org/10.1093/mnras/219.2.285.
- Walter Dehnen, *Phase-space mixing and the merging of cusps*, MNRAS 360 (2005) 892, https://doi.org/10.1111/j.1365-2966.2005.09099.x and https://arxiv.org/abs/astro-ph/0504246.
- Source static Bose halo: https://doi.org/10.1088/1475-7516/2021/11/055; our earlier papers 03,04,07,09.
