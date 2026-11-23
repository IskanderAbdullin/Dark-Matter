# Critical reversal: even NONoverlapping physical Bose DG cells can hit a false inverse-fugacity boundary, once boosted daughter angles are resolved (2026-11-23)

## New controlled falsification of the preceding partial repair

The [nonoverlapping physical p²dp entropy-dual energy-cell experiment](PHI-BOSE-CELL-ENTROPY-AUDIT.md) showed that if an **entire scalar cell has zero collision deposition**, its two local susceptibility coordinates have zero inverse-fugacity tangent. It passed an old failing implicit dt30 case on an N16 quadratic energy mesh. This was a scoped mathematical observation, **not** a proof that the discontinuous-P1 method always preserves Bose ηφ>0.

Independently, the [new exact boosted daughter log-Eφ angle](PHI-BOSE-DAUGHTER-LOG-ANGLE-AUDIT.md) showed that the historic ordinary nangle12 Gauss omitted ~10% of scalar inverse collision rate at pmax2 and severely under-sampled low-Eφ daughters. The next logical falsification was therefore mandatory: rerun the **physical moment tangent and actual nonlinear implicit time** with forward-resolved positive on-shell event weights and at least N16/N32 physical scalar energy cells. Retain a counterexample even if it contradicts an earlier apparent numerical repair.

## Exact algebraic no-go within EACH nonoverlapping cell

For two local scalar energy hats in a single physical positive-measure cell, the 2×2 physical Bose susceptibility J=dM/dη=−G is symmetric negative definite, G=[[a,b],[b,c]] with a,c>0, b>0, ac−b²>0. If both local scalar moments have **negative** collision source C=(−u,−v) with u,v≥0, then

    ηdot = J^{-1} C = G^{-1}(u,v),
    ηdot_left  = (c u − b v)/(a c − b²),
    ηdot_right = (a v − b u)/(a c − b²).

There is NO general sign constraint ηdot_left≥0. The explicit SPD numerical example a=3,b=1,c=2, u=.05,v=1 has C=(−.05,−1), yet **ηdot_left=−.18**, ηdot_right=+.59. Thus removing **cross-cell** overlap does not remove **within-cell** inverse susceptibility positivity obstructions. Neither the Pauli/Bose H theorem nor negative net scalar particle source proves the inverse fugacity lies inside its admissible tangent cone at the Bose boundary. A true condensate is a separate physical δ³(p) degree of freedom and is NOT implied.

## Actual physical Yukawa tree source: two scalar meshes, same well-resolved angle

Use the very same low-energy tree L=−yφ(barχ2χ1+h.c.) with on-shell χ1→χ2+φ and inverse Pauli/Bose, m1=.1 GeV,m2=.02 GeV,**mφ=86e−9 GeV**, y=.001, Tχ1=.9,Tχ2=Tφ=1 GeV, g1=4, physical parent pmax2, nrad24. For scalar ηφ(E), use positive quadratic E cells with TWO independent endpoints per cell, physical p²dp integrated moments; χ1/χ2 energy cells uniform. In each microreaction Γrest m1/E1 and boosted E2+Eφ=E1 exactly, full F/R = f1(1−f2)(1+bφ)−(1−f1)f2 bφ. The CM angle is resolved via the exact y=log(Eφ/Emin) transform, not the original low-order uniform cosine.

At nangle48, the *same* scalar physical weak source has Rφ≈−1.2954539554×10⁻¹¹ and microscopic H≈+1.8693373242×10⁻¹², independent of N16/N32 thermal affine mode reconstruction to floating precision. Every φ discrete source component is **nonpositive**, charges and energy exact to roundoff. Nonetheless:

| φ quadratic energy cells | first-cell source (left,right) | dηφ(rest)/dt | local ηrest/abs(dηrest/dt) | actual fully implicit dt30 |
|---:|---|---:|---:|---|
| N16 | (−1.335118e−12, −1.398127e−12) | **+5.837013e−9** | no rest-node bound | **passes** |
| **N32** | (**−6.574561e−14, −3.234469e−13**) | **−3.522727e−9** | **≈24.41291** | **fails closed** |

At N32 the actual implicit solver **accepts dt20** with min ηφ≈1.55523×10⁻⁸, positive finite physical entropy increase **+3.73863163×10⁻¹¹**, true Bose state, dark charge/energy residual around 4×10⁻¹⁵. It **rejects dt30 and dt1000** without flooring/clipping ηφ. The predicted linear local boundary dt≈24.41 is a diagnostic, not a proof that no other admissible root exists. A separate source-integrated physical moment finite-difference control reproduces J^{-1}C to absolute error <10⁻¹⁶ for a short tangent.

**Independent angular grids confirm the negative isn't merely a low-cosine artefact:** at the same N32, mapped nangle12/24/48/96 have dηrest/dt≈−3.4404/−2.9785/−3.5227/−3.1878 ×10⁻⁹, with linear limits ≈25.00/28.87/24.41/26.98, while the total physical Rφ and H remain angular-converged to much tighter relative precision. The **mode/first-cell** tangent is not uniformly stable under the same angular quadrature refinement.

The earlier ordinary-cosine nangle12 N32 test, by contrast, has **first-cell C=[0,0]** and **ηdotfirst=[0,0]** and admits dt30/1000. It also has an undercounted source Rφ≈−1.16493×10⁻¹¹, versus mapped −1.29545×10⁻¹¹. Thus the prior apparent first-cell safety partly relied on a **missing physical low-energy daughter angle region**; the fully mapped weak source exposes the nonoverlap method's **new and independent within-cell** obstruction.

## Nonmonotone mesh dependence is a genuine retained negative

With parent cutoff2, mapped nrad24,nangle48 and exactly the same Fermi/Bose **analytic thermal** state, dηrest/dt is positive at N8/12/16/20, negative at N24/32/40 and becomes precisely zero at N48/64 and finer because the quadrature events no longer touch the much narrower first cell. This is strong evidence of a grid/event-intersection pathology at extremely light mφ; it is **not** a converged statement about the continuous quantum collision at p=0. In particular, angular convergence of the *integrated* weak φ rate is much easier than convergence of the near-critical **mode-resolved** collision tangent.

A finite-time H theorem (N/E/S) and a positive Jacobian in physical integrated moments **do not imply** every convex finite-dimensional Bose manifold supports all physical source tangent directions. Simply changing the spatial support of hats or adding a post-hoc Bose η floor is not a valid proof of physical Bose production. A true zero-mode equation, if physics requires it, must be separately derived and source-normalized; otherwise a numerical false BEC will reappear as an artifact.

## Reproduction, research status and specific next hypothesis

- [Executable physical N16/N32/dt/angular counterexample](phi_bose_dg_tangent_counterexample.py).
- [Six focused regression tests](../tests/test_phi_bose_dg_tangent_counterexample.py).
- [Machine data](results/phi-bose-dg-tangent-counterexample-2026-11-23.json).
- [Dated stage](../archive/stages/2026/11/2026-11-23-phi-bose-dg-tangent-negative.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_dg_tangent_counterexample.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_dg_tangent_counterexample

Six focused tests passed locally in 13.03 s after avoiding redundant costly failed Newton attempts. The old DG success and angular correction code are separate prior commits; their CI status is tracked independently, not assumed. Published methodological context: [entropy-stable DG Boltzmann moments](https://arxiv.org/abs/1602.01312) and [quantum Boltzmann–Nordheim numerical BEC analysis](https://www.sciencedirect.com/science/article/pii/S0021999123002929), **neither** proves this new Yukawa φ scheme.

**Next scientific gate:** investigate physically admissible non-overlap **cell-average Bose entropy-dual with one scalar DOF per truly disjoint cell** or a positivity-constrained convex moment inversion, so a cell whose net source removes particles cannot require a negative rest-node η tangent from intra-cell interpolation. But a one-DOF constant η per energy cell does not reproduce energy exactly: a strictly energy-conserving paired source needs an additional independent reconstruction/constraint or dynamic cell-energy moment, or a separate higher-order entropy-compatible mechanism. Explicitly quantify that tradeoff before accepting a candidate. Resolve parent p/c/UV and nonthermal IR mode/time/zero-boundary at 3+ resolutions and only then integrate full χ1χ2–lepton vector Pauli, e±γν/Zprime finite plasma and FRW. **No physical Bose condensate, 86-eV relic ΩDM, Γelastic, frozen f(q), or CLASS-ready status.**
