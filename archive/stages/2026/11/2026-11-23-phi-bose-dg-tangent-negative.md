# Stage 2026-11-23 — nonoverlapping Bose DG still has a within-cell positivity obstruction

- **ID:** DM-KINETICS-20261123-BOSE-DG-WITHIN-CELL-TANGENT-NEGATIVE.
- **UTC launch beginning:** 2026-11-23T06:50:13Z. Exact stage finish can be reconstructed from Git commit; no invented durations.
- **Starting branch/HEAD:** research/kinetic-equilibrium-audit-20261121, [6a0d40a](https://github.com/IskanderAbdullin/Dark-Matter/commit/6a0d40a851c5b1e1d2a961b2d2a56c965933e12d), draft PR #1; preceding [real Scientific CI #45128900651](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45128900651) passed on the boosted daughter-log-angle stage, not on this new negative experiment.
- **Project context:** [AGENTS](../../../../AGENTS.md), [GOAL](../../../../GOAL.md), [PLAN](../../../../PLAN.md), [kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md), [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md), [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Question and independent falsification

Prior nonoverlap physical p²dp P1 cellwise η_Bose removed a **cross-cell** Gram-induced rest-node drift if there were no events in the first cell, and the N16 quadratic experiment passed previously failing dt30. But exact boosted Bose daughter angular integration recovered a low-Eφ source missing in the original underresolved uniform-c Gauss. Is the numerical repair universal? Test exact positive-η source, physical Γrest, same tree quantum Pauli/Bose collision, matched angular n12/24/48/96, scalar bins16/32/48/64 and actual fully implicit dt20/30/1000.

## Proof of within-cell inverse-fugacity failure

Even in a strictly **disjoint** cell, the two local linear endpoint hats overlap within the cell. With its SPD positive physical Bose Gram G=[[a,b],[b,c]], b>0, J=dM/dη=−G and both **negative** local source moments C=(−u,−v), u,v≥0, the local inverse-fugacity velocity is ηdot=G^{-1}(u,v), so ηdot_left=(c*u−b*v)/(a*c−b²) may be strictly **negative**. An exact counterexample G=[[3,1],[1,2]], C=(−.05,−1) gives ηdot_left=−.18. Neither negative net Bose scalar source nor H≥0 guarantees ηdot_rest≥0. It is a **numerical feasibility** obstruction, not condensation.

## Physical on-shell experiment and results

Code [kinetics/phi_bose_dg_tangent_counterexample.py](../../../../kinetics/phi_bose_dg_tangent_counterexample.py), [6 independent focused tests](../../../../tests/test_phi_bose_dg_tangent_counterexample.py), [full report](../../../../kinetics/PHI-BOSE-DG-TANGENT-COUNTEREXAMPLE-AUDIT.md), [compact data](../../../../kinetics/results/phi-bose-dg-tangent-counterexample-2026-11-23.json).

Inputs mχ1=.1 GeV,mχ2=.02GeV,mφ=86e−9GeV=86eV, Yukawa y=.001, thermal initial Tχ1=.9 GeV,Tχ2=Tφ=1GeV, parent pmax2GeV,nrad24, **mapped logarithmic daughter** angle nangle48, quadratic physical scalar energy cells, same local hats in Pauli/Bose forward/reverse interpolation and deposition, integrated p²dp Bose moments, positive event weights. The net weak φ rate Rφ=−1.295453955424128e−11 and microscopic H=+1.869337324222073e−12; all scalar moment source components are **nonpositive**; N/E per event conserved.

| Scalar cells | φ first-cell C(L,R) | ηdot_rest | actual implicit dt30 |
|---:|---|---:|---|
| N16 | (−1.335118e−12, −1.398127e−12) | +5.837013e−9 | passes |
| **N32** | (−6.574561e−14, −3.234469e−13) | **−3.522727e−9** | **fails closed** |

At N32, ηrest initially 8.6e−8, local linear boundary ~24.4129 GeV−1; **dt20 passes** with minimum ηφ=1.55523e−8, ΔS≈+3.739e−11, charge/E residual ~4e−15, but **dt30 and dt1000 reject** via physical positive-η Newton. This local linear bound is a diagnostic, not a proof that no other mathematical admissible root exists. Independently nangle12/24/48/96 give N32 negative ηrest_dot around −3.4404/−2.9785/−3.5227/−3.1878 ×10−9, despite nearly converged integrated Rφ/H. Old **uniform** cosine nangle12 had Cfirst=[0,0], ηdotfirst=0 and dt30 passing, but its integrated φ rate was ≈10.08% too small from missing boosted forward daughter events.

On nangle48 the **rest η tangent changes sign across bins**: N8/12/16/20 positive, N24/32/40 negative, N48/64 and finer zero because the finite event quadrature misses a narrowing first cell. This nonmonotone mode-grid behavior is not physical continuum proof. All six focused tests passed in a clean local sandbox (~13 seconds), no hidden clipping or invented scalar zero mode. New Scientific CI to verify separately after commit.

## Scientific decision and next executable stage

**Negative result preserved:** previous DG nonoverlap partial solution solves a zero-first-cell cross-cell coupling problem **only if no events**; it does **not** solve the general within-cell positive Bose manifold tangent. The physical weak rate can be angular-converged while first-cell source is not. Do not infer BEC from failures or repair by flooring η or inventing δ(p).

Next P0: compare a truly one-parameter nonoverlapping **cell average** entropy-dual Bose state (source sign then determines η direction) against an exactly energy-conserving model with additional independent energy variable, or a physically constrained moment inversion; derive and test quantum detailed balance, number/E/H before using in time. On ≥3 q/parent-p/angular/UV/dt resolutions, preserve full negative source and any spectral loss. Full χ1χ2+SM/eγν/Zprime/FRW and absolute Ω/BEC/frozen f(q)/CLASS not established.
