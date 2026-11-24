# Negative: quadratic Bernstein Bose moment projection is NOT pointwise sign-preserving (2026-11-24 UTC)

## Scientific question

The physical χ1(.1 GeV)↔χ2(.02 GeV)+φ(86 eV) Yukawa source at vacuum y=.001, parent FD T1=.9 GeV, daughter FD and scalar BE T=1 GeV, all chemical μ=0, **depletes each accessible scalar mode** for every finite parent P. This follows independently from detailed balance:

    F/R = exp[E1(1/Tdaughter−1/Tparent)] = exp(−E1/9) < 1

for every actual event E1>0, so pointwise dbφ(k)/dt<0 and true dηφ(k)/dt=−b_dot/[b(1+b)]>0 in open on-shell support. The previous DG-P1 affine entropy-variable moment scheme nonetheless could project g(k)=dotη(k)>0 into a negative endpoint derivative. Does a **nonnegative quadratic Bernstein-P2 basis** with an *exact physical support-edge split* fix this in general?

**No.** It greatly improves the tested P2 grids, but it does NOT guarantee pointwise positive projected dotη for every cell/cutoff. This is a scoped mathematical/numerical counterexample, not physical Bose creation, a negative phase occupation or a condensate.

## Independent physical experiment

Select the first physically **accessible** scalar energy bin [L,U] whose left edge is the exact massive on-shell kmin(P)=mφ cosh(max[0,asinh(kstar/mφ)−asinh(P/m1)]). The upper edge is the first original quadratic scalar grid energy U= mφ+[sqrt(P²+m1²)+.2−mφ]/N². Thus both the physically disconnected low-energy interval and the active interval are represented separately.

Use the nonnegative quadratic Bernstein hats B0=(1−x)², B1=2x(1−x), B2=x², x=(k−L)/(U−L). Independently integrate the **scalar-energy-first, two-sided-inverse-parent** physical modal source bdot(k)=S(k)−D(k)b(k) over physical DOS Dφ(k)=k sqrt(k²−mφ²)/(2π²). The three physical moments C_a=∫ Dφ B_a bdot and physical positive Bose susceptibility G_ab=∫Dφ B_a B_b b(1+b) give the unique projected entropy-variable tangent c=−G^-1 C. A second implementation with the common parent-angle **three-species physical event list** independently cross-checks the result. No synthetic mode capacities or guessed average energies enter.

The exact projected polynomial is g_h(x)=c0(1−x)²+2 c1x(1−x)+c2x². The exact interval minimum is checked analytically from its two endpoints and stationary point x*=(c0−c1)/(c0−2c1+c2), when this is a convex interior quadratic.

## Numerically converged counterexamples

The table uses independent physical scalar k-first × inverse-parent Gauss 80×80, with positive on-shell DOS. Values of c are in natural GeV inverse-fugacity rate units.

| Parent cutoff P (GeV) | Original Nφ | (c0,c1,c2) (GeV) | Minimum of PROJECTED g_h(k), GeV | Sign |
|---:|---:|---|---:|---|
| 2 | 8 | (+9.36137938e−9, **−1.41111862e−9**, +1.96315296e−9) | **+1.15832549e−9** | pointwise positive |
| 2 | 12 | (+1.02079265e−8,+3.88422219e−9,+2.34395549e−9) | +2.34395549e−9 | positive |
| **8** | **8** | (+1.29650241e−8, **−7.62967217e−9**, +3.23779803e−9) | **−5.15977428e−10** | **false negative** |
| 8 | 12 | (+2.00760739e−8,−8.74648969e−9,+4.51949766e−9) | +3.38160536e−10 | positive |
| 8 | 16 | (+2.52223095e−8,−6.43041690e−9,+5.13646100e−9) | +2.04081340e−9 | positive |
| **16** | **8** | (+8.39458629e−9,−5.58824193e−9,+2.22194588e−9) | **−5.77071638e−10** | **false negative** |
| **16** | **12** | (+...,−7.87230834e−9,+...) | **−4.31773673e−10** | **false negative** |
| 16 | 16 | (+1.89263642e−8,−8.49835385e−9,+4.26106014e−9) | +2.09643893e−10 | positive |
| 16 | 32 | (+2.69934485e−8,+6.44727095e−9,+4.46155531e−9) | +4.46155531e−9 | positive |

Most decisive: **P=8 GeV, Nφ=8**, the fitted quadratic has g_h,min=−5.15977428×10⁻¹⁰ GeV at x=0.654586, k=**0.08397890952 GeV** in the physically supported first cell. But the true microscopic scalar-mode source at the SAME k gives **g_true(k)=+9.68624932×10⁻¹⁰ GeV**. The physical susceptibility is well-conditioned (condition number ≈10.00065), and the negative g_h result persists as independent scalar Gauss and parent Gauss orders increase 24→32→48→64→80. It is not a Bose IR cutoff/no-support zero or floating-point noise.

**Important nuance:** c1<0 is *not* enough to establish a negative projected function: for P2 N8 and P8 N12/N16 the middle coefficient is negative but g_h(x)>0 throughout. In particular for c0,c2≥0, c1≥−sqrt(c0*c2) is a sufficient pointwise nonnegativity condition even if c1<0. Our negative P8 N8 obeys c1<−sqrt(c0 c2), and direct polynomial minimization proves an actual negative subinterval.

At finer Nφ16 (P8) and Nφ16/32 (P16), the shown first-supported-cell projected tangent returns positive, so mesh refinement *can* remedy the specified counterexample. It does not establish a universal global positivity bound across all cells/nonequilibrium states.

## Physical interpretation and reproduction

The real physical Bose spectrum retains positive occupation and true pointwise depletion; only a **finite-dimensional polynomial moment projection** gives the false sign. An entropy-stable implicit finite-time step can still remain admissible because a small negative local eta-dot does not force η(E) below zero; a global positivity theorem must not be inferred from good finite-time samples. A genuinely shape-preserving projected operator or an explicitly constrained moment cone may be required in more difficult regimes.

Source: [independent k-first and event-first tangent audit](phi_bose_bernstein_p2_tangent_counterexample.py). Tests: [five targeted physical regressions](../tests/test_phi_bose_bernstein_p2_tangent_counterexample.py). Commands:
- OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_bernstein_p2_tangent_counterexample.py
- OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_bernstein_p2_tangent_counterexample

Independent local Python/SciPy n=16/24/32/48/64/80 Gauss refinements show the sign and the printed numerical values; the new GitHub tests/source CI must be checked separately before publishing as independently verified project evidence.

**NEXT_ACTION:** audit and bound all-cell shape-preserving P2 energy projection for arbitrary nonthermal Bose initial spectra and cutoff P, and check the actual physically constrained finite-time implicit root. Physical IR/UV/time convergence, shared event nullspace, SM thermal masses/FRW and relic/observational likelihood remain unproved.
