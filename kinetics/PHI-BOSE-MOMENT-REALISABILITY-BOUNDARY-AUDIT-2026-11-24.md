# Exact physical Bose first-cell two-moment realizability boundary (2026-11-24 UTC)

## Question

Does the **physically converged weak** Yukawa reverse source really drive the positive-mass 86-eV scalar out of the space of nonnegative physical occupations, or only out of the **nonoverlapping DG-P1 affine inverse-fugacity two-moment ansatz**?

We answer this sharply for a *frozen collision-source moment ray*, not a dynamically recomputed strong collision operator. The prior source, its four weak N/Q moments, and the negative inverse-fugacity tangent have independent parent-shell on-shell convergence controls.

Physical benchmark: m1=.1, m2=.02, mφ=86e-9 GeV, Yukawa y=.001, g1=4, parent p∈[0,P], FD parent T1=.9, FD daughter and BE scalar T=1, zero chemical potentials; actual boosted χ1↔χ2+φ weak decay kinematics and exact scalar physical p²dp density of states. Scalar first cell is [mφ,U] on the existing quadratic energy grid at N=24,32,40 (with P=2 GeV unless stated).

## Rigorous finite two-moment face and uniqueness

Define two nonoverlapping local hats hL=(U−E)/(U−m), hR=(E−m)/(U−m), m=mφ, with hL+hR=1. For a smooth Bose spectral function b(E)>0, define physical number/energy-equivalent moments

M_a[b] = (1/2π²)∫_0^sqrt(U²−m²) p²dp h_a(E(p)) b(E(p)), a=L,R.

For the affine inverse fugacity family η(E)=ηL hL+ηR hR, b(E)=1/[exp(η(E))−1], admissibility requires ηL≥0, ηR>0; ηL=0 is the *critical face*, with finite p²dp number and energy for **m>0** and ηR>0 even though b(p→0) diverges ∝p^−2. (At p=0, p²b→2m(U−m)/ηR. This does **not** establish a condensate δ³(p).)

For interior ηL>0, ηR>0, the response G_ab=(1/2π²)∫p²dp h_a h_b b(1+b) is strictly positive definite. Since dM=−G dη, the map is locally one-to-one. At **fixed M_R**, varying ηL while adjusting ηR gives

dM_L/dηL |_(M_R fixed) = −(G_LL − G_LR²/G_RR) < 0.

Therefore, provided the boundary state at ηL=0 has the desired M_R, the maximum possible M_L in the **entire positive affine-η family** at that M_R is precisely M_L[ηL=0,ηR(M_R)]. A target with larger M_L cannot be attained by *any* positive affine inverse-fugacity parameters. This is an exact scoped feasibility inequality, not a mere solver failure.

Freeze the independently calculated signed physical weak sources in this first cell at the thermal initial state: C_L=(U Rφ−Qφ)/(U−m), C_R=(Qφ−m Rφ)/(U−m), both negative in the studied inverse state. Set M_target(t)=M(ηL=m/T,ηR=U/T)+t(C_L,C_R), with natural time t in GeV^−1. Find the **exact first critical face intersection** by solving two physical moments, M(0,ηR*)=M0+t* C, using adaptive p²dp quadrature over threshold-split physical momentum. The integral is evaluated as **M−M0 directly**, so large thermal M0 subtraction cannot hide the ∼1e−12 GeV³ effect.

## Independently reproduced result

| P (GeV) | first-cell N | ηL dot (GeV) | linear ηL boundary (GeV^−1) | exact frozen-source face t* (GeV^−1) |
|---:|---:|---:|---:|---:|
| 2 | 16 | +5.704290652e−9 | — | not approached locally |
| 2 | 20 | +2.710623342e−9 | — | not approached locally |
| 2 | 24 | −4.934728039e−10 | 174.27505 | **174.23628158** |
| 2 | 32 | −3.264074504e−9 | 26.34744 | **26.38466116** |
| 2 | 40 | −5.764160436e−10 | 149.19779 | **149.52832484** |
| 2 | 48 / 64 | 0 | — | physically forbidden first-cell source |
| 8, 16 | 16 / 24 / 32 / 40 | positive | — | not approached locally in these grids |

At N32,P=2 the **exact critical endpoint** ηL*=0, ηR*=0.0021512311069491253 (initial ηR0=U=0.0021509632982910924). The frozen physical sources are C_L=−6.581061116761629e−14 and C_R=−3.0943944628484195e−13 GeV⁴. The two separately computed moment ratios ΔM_L/C_L and ΔM_R/C_R agree to relative 5e−13. An independent 50-decimal mpmath p-integration check of the same physical boundary gives 26.38466116410699 and 26.38466116403818 GeV^−1 from the two coordinates: this independently supports the SciPy value, with differences at ∼1e−10 absolute caused by finite root precision.

**Rigorous exclusion after crossing for the frozen ray within this family:** at t=1.1t*, fixing the target M_R and choosing the boundary ηL=0 that maximizes achievable M_L gives

| N | M_L(target) − max allowed M_L (GeV³) at 1.1t* |
|---:|---:|
| 24 | +4.16867436e−13 |
| 32 | **+2.34632711e−13** |
| 40 | +1.50277307e−13 |

The M_R matching residual is <5e−22 GeV³. At t=0.8t* the corresponding gap is negative for all three N. This establishes a *real* crossing of the DG-P1 representable domain rather than a Newton or floating-point convergence issue. It is NOT a proof of finite-step failure when the physical collision operator is recomputed nonlinearly at every trial state.

## Constructive positive distribution with the SAME excluded moments

The frozen ray is **not outside the entire physical nonnegative Bose phase-space moment cone**. Define, with the *same* local hats,

b_new(E)=b_T(E)[1+a_L h_L(E)+a_R h_R(E)].

Use the independently physical DOS matrix A_ab=∫p²dp/(2π²) b_T(E) h_a h_b, positive definite. Solve A a=t C, so both physical moments **exactly** equal M0+t C, while b_new(E)>0 everywhere whenever both endpoint multipliers 1+aL and 1+aR are positive.

At t=2t* the explicit coefficients are:
- N24: (aL,aR)=(-0.00098339118,-0.00125584148);
- **N32: (+0.00015145189,−0.00032940754)**;
- N40: (+0.00035631035,−0.00039517961).

All endpoints are >0.9987, all target moments positive, and independent direct physical DOS moment reconstruction agrees within 5e−21 GeV³. Even at t=10t*, all three N admit such positive non-affine spectra; none need condensate injection or η-clipping. This **constructively distinguishes artificial affine-η moment-family failure from physical positive-distribution impossibility**. These non-affine trial spectra are not claimed to be the actual self-consistent Boltzmann time evolution or entropy-producing solution.

## Checks, method and scientific limits

- [Implementation](phi_bose_moment_realisability_boundary.py), [five focused tests](../tests/test_phi_bose_moment_realisability_boundary.py).
- Test command: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_moment_realisability_boundary.py
- Audit command: OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_moment_realisability_boundary
- Five local focused tests passed in 7.85 seconds, including exact 2-moment face, inside/outside signed gaps, explicit physically positive non-affine counterexample, parent P=8/16 sign controls, hard P=2 N48/64 support, stable E−m=p²/(E+m) threshold geometry. New public Scientific CI must be checked independently; test success certifies only these computations.

**Scope and caveats:** frozen prescribed thermal weak source, a *single* local scalar bin, a specific DG-P1 affine inverse fugacity family. No finite-step full χ1/χ2/φ self-consistent strong collision integration, no proof all physical spectra realize arbitrary finite-time sources, no physical Γ(T)/SM heat bath/FRW, no true 86-eV zero mode or condensate, no Ωφ, CLASS or observational likelihood. N=48/64 first-bin zero arises from the hard P=2 shell support floor, not Bose saturation. A critical endpoint ηL=0 with integrable Bose infrared number is not in itself a dynamically formed condensate.

**NEXT_ACTION:** move beyond affine-η two-DOF family: propose a physical positive multi-moment spectral representation (or constrained entropy dual with additional local degrees) that matches independently converged on-shell N/Q, exact two charges/three Q event ledger and H-theorem; independently test full nonlinear strong F/R at each microevent and finite-time positivity without clipping. Audit true graph nullspace, IR/UV/time convergence and finite SM+FRW before physical relic claims.
