# Physical-shell resolved Bose DG-P1 tangent cone: robust negative result (2026-11-24)

## Physical question

Earlier nonoverlapping Bose DG-P1 calculations found a **negative first-node inverse-fugacity tangent despite negative Bose particle source**. Those calculations used less well-converged boosted scalar angles. The new independent exact shell-crossing physical Yukawa weak integration can now isolate whether the negative tangent is a true property of this *numerical moment manifold* or a sparse-angle artifact. Here the parent-cutoff is P=2 GeV; m1=.1, m2=.02, mφ=86e−9 GeV, Yukawa y=.001, g1=4, parent FD T=.9 and daughters FD/BE T=1, all chemical potentials zero.

This is an **instantaneous tangent**, NOT a finite-step nonlinear existence/uniqueness or a physical condensate calculation.

## Derivation using only independently converged first-bin R and Q

For a physical scalar energy bin [L,U], choose two local entropy-dual affine hats hL=(U−E)/(U−L) and hR=(E−L)/(U−L), with no cross-bin supports. The exact physically integrated weak moment sources can be recovered *solely* from the true on-shell bin number and scalar energy sources:

\[
C_L=\frac{U R_\phi-Q_\phi}{U-L},
\qquad
C_R=\frac{Q_\phi-LR_\phi}{U-L}.
\]

No new daughter-angle sampling or guessed bin-center occupations enter this projection. For Bose inverse fugacity eta(E)=E/T (no condensate), independently integrate the physical positive susceptibility

\[
G_{ab}=\frac{1}{2\pi^2}\int_{p(L)}^{p(U)}
p^2\,dp\, h_a(E)h_b(E)\,b(E)[1+b(E)].
\]

The Jacobian dM/deta=−G is negative definite, so the **unique instantaneous tangent** obeys \(\dot\eta=-G^{-1}C\). The 2×2 matrix is evaluated by adaptive SciPy quadrature over physical p, deliberately split near p≈mφ to resolve the extremely light massive threshold, independently of the parent/daughter reaction quadrature.

## Numerical results

All results use exact parent shell-crossing splitting and positive log-daughter-angular quadrature, from the preceding physical N/Q benchmark; the reference 10 parent nodes *per smooth segment* ×16 log-angle nodes has checked convergence against 6×8,16×32,32×48.

| Number N of quadratic scalar bins at P2 | physical first-bin Rφ (GeV⁴) | dη_left/dt (GeV) | local linear boundary time (GeV⁻¹) |
|---:|---:|---:|---:|
| 8 | −5.560284984612e−12 | +6.301800104e−9 | none from first node |
| 12 | −3.912375560091e−12 | +7.201265406e−9 | none |
| 16 | −2.762149360447e−12 | +5.704290652e−9 | none |
| 20 | −1.898871853310e−12 | +2.710623342e−9 | none |
| 24 | −1.236651580795e−12 | **−4.934728039e−10** | 174.275 |
| 32 | −3.752500574525e−13 | **−3.264074504e−9** | **26.3474** |
| 40 | −2.320689918353e−14 | **−5.764160436e−10** | 149.198 |
| 48 | exactly 0 | 0 | no physically accessible first-bin event |
| 64 | exactly 0 | 0 | no physically accessible first-bin event |

At N32, both true scalar moment sources are **negative**:
CL=−6.58106111676e−14, CR=−3.09439446285e−13 GeV⁴, yet the first inverse-fugacity derivative is negative. This sign is not numerical noise: 6×8,10×16,16×32,32×48 shell-partitioned angular/radial grids agree on dηleft to relative ≪1e−9 (limited by floating point and fixed physical susceptibility). Exact positive G eigenvalues independently certify the invertible local map. For N32, its physical matrix is approximately

G = [[3.63190856e−5, 1.81607833e−5],
     [1.81607833e−5, 3.63215638e−5]].

An independent high-temperature near-massless asymptote gives \(G\simeq [U/(2\pi^2)] [[1/3,1/6],[1/6,1/3]]\) for mφ≪U≪T; numerical G differs from this simpler integral by less than 1.5e−4 relative. This controls the physical susceptibility separately from the collision source.

The earlier angular n48 N32 result dηleft≈−3.522727e−9 had about a 7.34% larger magnitude than this **parent-partitioned** value. The qualitative negative survives improved quadrature; its previous numerical coefficient was not fully converged.

## Interpretation and negative control

At N32 the thermal rest-node eta≈mφ/T=8.6e−8, giving a **linearized** tangent limit eta/|deta/dt|≈26.35 GeV⁻¹. This is a *local diagnostic* of the DG-P1 Bose positivity manifold; **not** a proof that every implicit finite-step solver fails above dt=26.35 or that a physical delta-condensate nucleates. The underlying source removes φ from every tested first-bin hat; it never requires physical Bose creation to generate this numerical obstruction. This refutes any unqualified claim that exact on-shell energy, weak entropy and negative Bose source guarantee a positive affine inverse fugacity tangent.

The N-dependence is highly nonmonotone and at N48/64 the first cell becomes **physically inaccessible at finite P=2**, not a converged continuous pφ=0 statement. Refining the scalar grid alone can therefore make the first-bin source disappear for purely kinematic reasons. The finite parent cutoff cannot access the true 86-eV zero-momentum region, and a δ³(p) condensate needs an independent derived kinetic channel.

## Reproduction

- [Standalone tangent and physical susceptibility code](phi_bose_exact_cutcell_tangent.py)
- [Four tests](../tests/test_phi_bose_exact_cutcell_tangent.py)
- [Underlying independently validated positive physical weak cut-cell N/Q](PHI-BOSE-PHYSICAL-CUTCELL-NQ-AUDIT-2026-11-24.md)
- Commands: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_exact_cutcell_tangent.py
- Commands: OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_exact_cutcell_tangent

Four local tests passed in 6.82 seconds, including the sign and convergence over 4 quadratures for N16/24/32/40, true N48/64 hard support zero, exact Jdeta=C, and independent high-T asymptote of physical susceptibility. Public Scientific CI of this new stage must be verified separately.

**NEXT_ACTION:** now that the negative tangent is accurately established, test an explicitly *realizable and entropy-stable* physical two-moment finite-cell step (e.g. convex constrained moment update or exact local quantum microevent extents with shared phase capacity). Before accepting it require physical weak N+three-Q agreement, detailed balance, Pauli/Bose constraints, both charges, total energy, time/mesh/IR tests. Do not floor eta, post-inject condensate or rescale Ω to hide a failure. Strong self-consistent χ+SM+FRW relic/CLASS/observations remain unestablished.
