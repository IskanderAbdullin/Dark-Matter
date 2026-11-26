# A Bose IR dual-moment tangent-cone obstruction despite strictly negative physical scalar collision sources (2026-11-22)

## Question and falsification of an unsafe BEC inference

A physically integrated p²dp inverse-fugacity χ1↔χ2+real φ Galerkin solver gives exact energy, two dark charges and a quantum H theorem. Nonuniform E-hats dramatically improve the source of a sharply curved nonthermal Bose distribution, but at mφ=86eV and a thermal massless-limit-like initial ηφ(mφ)=mφ/T≈8.6×10⁻⁸, a positive damped Newton implicit time step fails on some E grids. Could this be a real condensate, or a **purely numerical incompatibility between a projected collision moment direction and the admissible ηφ≥0 manifold**?

The physical Yukawa event source is evaluated with a cold χ1 (FD Tχ1=.9), warm χ2/φ (Tχ2=Tφ=1), chosen but genuine tree Γrest (mχ1=.1GeV,mχ2=.02GeV,mφ=86e−9GeV,y=.001), boosted on-shell 1↔2 and full Pauli/Bose F−R. Every physical φ event net flux is **negative** under these prescribed thermal occupations. The discrete φ moments in the affected first IR hat are **exactly zero source**, while all other φ source components are nonpositive and global event entropy positive. No Bose zero-mode δ³(p) is introduced.

## Exact two-overlapping-hat algebraic counterexample

Let two positive overlapping piecewise-linear Bose energy hats have physical moment susceptibility

    J = ∂M/∂η = −G,    G = [[a,b],[b,c]],    a,c>0, b>0,
    D = ac−b² > 0.

The positivity and overlap of hats guarantee symmetric negative definite J, with b>0. Suppose the legitimate projected collision source removes population from the second hat but not the first:

    C = (0,−q),   q>0.

The only instantaneous inverse-fugacity velocity satisfying dM/dt=C is

    dη/dt = J⁻¹ C = (−bq/D, +aq/D).

**The first Bose inverse fugacity DECREASES although its own moment source is zero and the only physical event removes bosons from the other hat.** This occurs because overlapping energy hats have a non-diagonal physical mass/susceptibility matrix: to preserve M_first while lowering its shared neighbor, the reconstructed first η must move downward. At a critical boundary η_first=0 this velocity points OUTSIDE the admissible positive-inverse-fugacity tangent cone. At a small positive η0, the first-order estimated time to the boundary is

    t_limit,lin = η0 / (bq/D).

This is an exact mathematical counterexample to the claim that quantum positivity of event rates and the H theorem automatically guarantee **positive inverse-fugacity nodal time evolution** for continuous overlapping Galerkin hats. It does **not** imply a physical Bose condensation phase transition; all collision events can be depleting φ.

Independent tiny arithmetic control G=[[3,1],[1,2]], q=1 yields dη/dt=(-.2,+.6), exactly matched by a separate NumPy solve.

## Physical finite Yukawa event/source/Jacobian measurement

Compute Jφ from the physical p²dp integral of the nonlinear Bose η(E)=Σh_iη_i, and Cφ from the SAME physical vacuum-width, moving-parent Pauli/Bose events (mother pmax2, nr12,nangle12). Independently solve Jφ dηφ/dt=Cφ at N8/16/32/64 and four uniform/quadratic/cubic/two-zone η-hat grids, tracking total dark charges, E, event H, all per-hat scalar net source signs and the smallest linearized η>0 crossing time.

| φ grid | dηφ,rest/dt | First scalar hat source | Linearized earliest Bose-boundary time |
|---|---:|---:|---:|
| N8 uniform | +1.444e−9 | negative | ~3.20e9 (another bin) |
| N8 cubic | **−2.327e−9** | 0 | **36.96** |
| N16 uniform | +2.198e−9 | negative | ~6.36e10 (another bin) |
| **N16 quadratic** | **−2.998203e−9** | **0** | **28.68385** |
| N32 cubic | **−3.915e−11** | 0 | **2196.81** |
| N32 two-zone | **−1.188e−9** | negative | **72.39** |
| N64 cubic | +8.15e−14 | 0 | ~7.30e6 (another bin) |
| N64 two-zone | **−1.263e−9** | 0 | **68.11** |

All 16 true physical on-shell source calculations have negative net φ rate and positive total microscopic entropy H≈+1.64846911e−12; no +φ first-hat source is needed to drive ηrest negative in the inverse moment representation. Total charge and energy source errors remain roundoff-level. The physical scalar mass is not itself a measured DM property — 86 eV is a historical numerical benchmark.

The N16 quadratic sensitivity with ηφ,rest≈8.6e−8 predicts a **linear** admissibility duration ≈28.68 in synthetic GeV⁻¹. The actual independently solved implicit moment update **succeeds at dt=20** with min ηφ≈2.60365e−8 and physical entropy increase, but **fails closed at dt=30** because damped Newton cannot find an improving residual while keeping all ηφ>0. That correspondence diagnoses a real discretized tangent issue, not merely a random numerical tolerance. A source-level finite difference test verifies Jφ⁻¹Cφ recovers the physical moment increment to absolute <1e−16 for an admissible small step.

Changing only the hypothetical φ rest mass in the SAME N16 quadratic test from 86eV to 8600eV increases the predicted linear bound from 28.684 to 2866.44 (≈100×), while first-node sensitivity dη0/dt remains about −3e−9. Further mφ=860eV gives ~286.90. Thus the numerical obstruction gets more severe as a very light mass approaches the ηrest=0 boundary, exactly where the former false-condensate rest-node collocation was already pathological.

## Consequences, negative controls, next candidate

**Rigorous conclusion within this numerical formulation:** a physical positivity/entropy-compatible microscopic χ1↔χ2+φ collision source can generate a semi-discrete Galerkin moment derivative outside the positive inverse-fugacity tangent cone at a Bose rest-node boundary. Strict H, dark charges, on-shell energy and correct integrated equilibrium number do not rule it out. For positive ηrest the linear bound is only a local predictor, not a proof that every implicit method must fail for dt above that value. The actual failures at dt≥30 are separately observed for this solver.

**Not established:** existence of a true dynamically generated BEC, physical zero-mode nucleation, absolute 86-eV relic density, realistic e±γν/Z′/FRW kinetic rates, or an infrared-complete χ1↔χ2+φ mode-dependent time solution. Do not solve this by silently flooring φ η or inventing δ³(p) number.

An explicit structural alternative to investigate is **non-overlapping discontinuous energy-cell entropy-dual moments** (within each cell two local energy hats for N and U, no across-cell sharing). If the first physical cell has no events, both local moment sources vanish, and its two η parameters can remain unchanged, potentially removing this particular inverse-mass-matrix cross-cell sign. This idea requires an independent on-shell F/R/H proof and actual positive implicit time validation: discontinuity at boundaries, the cost of extra modes, and physical q resolution might introduce other problems. Also test independent pmax12/16 and 3+ p/c/IR/dt resolutions before treating any result as a cosmological rate.

## Reproduction and versioning

[Exact two-hat formula and physical 16-grid tangent audit](phi_bose_moment_tangent_gate.py); [four focused tests](../tests/test_phi_bose_moment_tangent_gate.py); [primary 16-grid/4-mass JSON](results/phi-bose-moment-tangent-gate-2026-11-22.json); [permanent stage](../archive/stages/2026/11/2026-11-22-phi-bose-moment-tangent.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_moment_tangent_gate.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_moment_tangent_gate

Four focused tests pass locally. Real GitHub Scientific CI should be checked independently after the commit. The incoming Yukawa Γrest is conditionally tree-physical but masses/T/initial occupations are illustrative, and no actual cosmological ΩDM, BEC, Γelastic, frozen f(q) or CLASS-ready result is established.
