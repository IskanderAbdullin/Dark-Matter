# Stage 2026-11-22 — fully implicit, positive Bose/Fermi thermodynamic Yukawa moment step

**ID:** DM-KINETICS-20261122-PHI-BOSE-IMPLICIT-MOMENT-STEP. **Session UTC:** started 2026-11-22T21:20:21Z from isolated sandbox clock, exact substage elapsed not fabricated. Existing public Dark-Matter research branch/draft PR #1, immediately after [physical scalar Bose moment stage](2026-11-22-phi-bose-moment-galerkin.md), not a duplicate project. [GOAL](../../../../GOAL.md), [PLAN](../../../../PLAN.md), [kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md), [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md), [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Falsifiable question, method, sources

New physical p²dp moment basis has exact quantum instantaneous on-shell χ1↔χ2+φ conservation/H and correct scalar number (earlier 86eV collocation false ×87199/1363/22). Challenge: **finite-step evolution** without losing Bose positivity, Fermi bounds, N/E or H. Derived fully implicit backward Euler on integrated moments M(ηnew)−M(ηold)=dt C(ηnew), analytic symmetric-negative susceptibility and analytic event Jacobian, safeguarded Newton line search ηφ>0 with no clipping/reset or post-step energy repair. Concavity of true thermodynamic entropy proves Snew−Sold≥ηnew·ΔM=dt∑W(F−R)ln(F/R)≥0. An actual positive root and residual are required for the claim; bφ=0 and true condensate δ³(p) are not covered.

Files: [detailed proof and numerical negative](../../../../kinetics/PHI-BOSE-IMPLICIT-MOMENT-TIME-AUDIT.md), [solver](../../../../kinetics/phi_bose_implicit_moment_step.py), [9 focused tests](../../../../tests/test_phi_bose_implicit_moment_step.py), [compact numbers](../../../../kinetics/results/phi-bose-implicit-moment-step-2026-11-22.json). Run: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_implicit_moment_step.py tests/test_phi_bose_moment_galerkin.py ; OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_implicit_moment_step --audit. Nine new + eight prior tests passed locally. New repo Scientific CI to verify after commit; no broad unnecessary CI or private Actions requested.

### Positive and retained negative results

At mφ=86e−9 GeV, 16 energy bins and 16 synthetic on-shell positive-weight events, implicit dt=.25/.5/1/2 uses 4/4/4/5 Newton iterations, physical ΔS=+1.28/+2.00/+2.79/+3.48×10⁻⁶ exceeding dt·H, total N/E and charges at floating precision, Bose ηφ>0. First-order temporal self-convergence at t2, 2/4/8/16/32/64 steps, final adjacent-halving ratio≈1.964. Chemical fixed point stationary. Jacobian independently finite-difference checked.

**Critical negative:** with those sparse events, t2 and 8 implicit steps at N8/16/32/64 produce ΔNφ≈−5.993/−6.262/−2.830/−1.645×10⁻⁵: strong mode-grid nonconvergence even though entropy/charges/thermal initial N correct.

First independent physical matrix-element rate gate: for **chosen** scalar Yukawa y=.001, mχ1=.1,mχ2=.02,mφ=86e−9 GeV, exact tree Γrest=2.7501974166×10⁻⁹ GeV; on-shell p²dp Γrest m1/E1 dc/2 at pmax2, nr12,nangle12 yields Rφ=−1.16493012863×10⁻¹¹, positive H, matching a separate moving Yukawa weak code to better than 2×10⁻¹² relative. Physically weighted t=1e7 GeV⁻¹ nonlinear N8/16/32 ΔNφ≈−8.52569/−8.06803/−7.77340×10⁻⁵, still ≈9.7% N8→32 discrepancy: **rate input matching ≠ nonthermal strong convergence**. Synthetic temperatures, coupling and g-counting are not experimentally constrained.

## Scope, blockers and next executable scientific gate

No χ/SM lepton elastic collision in the **same** new moment solver, finite plasma/mediator/thermal masses, FRW H, bφ=0 critical/δ-mode, physical ΩDM/BEC/frozen f(q)/CLASS-ready. No physical Γelastic. P0 next: replace sparse mother event samples by joint 3+ incoming p/c/IR/mode/dt convergence with physical Γrest and strict entropy; then couple existing two-χ shared-lepton Pauli strong to identical moment/entropy formulation and real SM/Z′/FRW. Keep prior first-node pseudo-condensate evidence as permanent negative.
