# 2026-11-23 — Nonoverlapping physical Bose/Fermi entropy-dual cells: scoped IR realisability repair

- **ID:** DM-KINETICS-20261123-DISCONTINUOUS-P1-BOSE-TANGENT
- **Active launch started:** 2026-11-23T06:50:13Z UTC, recorded directly by isolated sandbox clock; final elapsed time must be measured, not inferred.
- **Starting HEAD:** [01e93ff](https://github.com/IskanderAbdullin/Dark-Matter/commit/01e93ff984021e4514d1c7564d8dff3419c10b57) in public research/kinetic-equilibrium-audit-20261121, draft PR #1. [Actual Scientific #45074692117](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45074692117) passed for the **preceding** overlapping-hat tangent negative study; new-stage CI separately pending.
- **Task context:** [AGENTS](../../../../AGENTS.md) → [GOAL](../../../../GOAL.md) → [PLAN](../../../../PLAN.md) → [kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md) → [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md) → [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Question and predeclared criteria

The preceding [IR moment-tangent negative stage](2026-11-22-phi-bose-moment-tangent.md) showed a **false Bose boundary tangent obstruction** from the overlap of energy hats: even with zero first-hat φ collision deposition and all scalar rates negative, the inverse susceptibility couples adjacent nodes so ηφ(rest) drops toward zero, causing the old implicit moment solver to reject dt30 on N16 quadratic grid (dt20 passed). Test the **structural** alternative: discontinuous cellwise linear energy basis, two independent endpoint Fermi/Bose inverse fugacity coordinates per physical energy cell.

Criteria: physical p²dp integrated BE on 3 meshes; exact N/E, Pauli/Bose entropy with identical local hats for reconstruction/source; no cross-cell susceptibility; common-T nonzero chemical fixed point; no artificial first-cell drift with zero local rates; nonlinear implicit dt20/dt30/dt1000 without clipping. Preserve any q-grid drift and negative results.

## Derivation and method

Per energy cell [Ej,Ej+1], use local hats a=(Ej+1−E)/ΔE, b=(E−Ej)/ΔE with zero support outside this cell and duplicated endpoints at interfaces. Physical moment Mja=g/(2π²)∫cell p²dp hja(E(p)) n(E(p)), U=Σedge Eedge M. Fully boosted positive-weight Yukawa Pauli/Bose event is evaluated/deposited with **the same hats**, so N1+N2 and N1+Nφ and total energy conserve; η·C=ΣW(F−R)ln(F/R)≥0. The scalar b=1/(expη−1) requires strict ηφ>0. Actual Bose N/S integrated in momentum, not by evaluating BE(mφ) at a finite-volume rest node.

Each physical per-cell susceptibility J=dM/dη=−∫cell dΠ h_a h_b n(1±n) is negative-definite and has exactly zero cross-cell blocks. If an entire first cell is untouched, Cfirst_left=Cfirst_right=0 and **both** first-cell inverse fugacity tangent derivatives are zero. Previous overlapping continuous hats did not guarantee this. Fully implicit Mnew−Mold=dt Cnew with positive-η Newton/backtracking and measured finite entropy is used, with no post-step clipping.

## Source, tests and measured evidence

- [Self-contained implementation](../../../../kinetics/phi_bose_cell_entropy.py).
- [Seven focused tests](../../../../tests/test_phi_bose_cell_entropy.py).
- [Detailed scientific audit](../../../../kinetics/PHI-BOSE-CELL-ENTROPY-AUDIT.md).
- [Four-grid/time JSON](../../../../kinetics/results/phi-bose-cell-entropy-2026-11-23.json).

Reproduce:
    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_cell_entropy.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_cell_entropy

Seven tests passed **locally** in 0.92 seconds with NumPy/SciPy; GH Scientific CI for stage remains to be confirmed. Toy parameters: m1=.1 GeV,m2=.02 GeV,mφ=86e−9 GeV,y=.001,T1=.9 GeV,T2=Tφ=1 GeV, on-shell positive tree Yukawa parent pmax2,nr12×nangle12,g1=4. Scalar quadratic energy-cell families N8/16/32/64.

At N16 the independent φ first-cell collisional deposition [0,0] gives inverse fugacity derivatives [0,0]. Old overlapping first-hat C=0 nevertheless gave ηrest_dot≈−2.998×10⁻⁹ and implicit dt30 failure. New finite-time N16 dt20/30/1000 all pass, ηφ(rest)=8.6×10⁻⁸ unchanged, dark charges/energy up to ~10⁻¹⁴. At dt30 N16 ΔNφ=−3.494750×10⁻¹⁰ and physical ΔS=+4.945377×10⁻¹¹; dt1000 ΔNφ=−1.164681×10⁻⁸ and ΔS=+1.648143×10⁻⁹. Initial source Rφ=−1.16493012863×10⁻¹¹, H=+1.64846911074×10⁻¹² agrees with prior moving Yukawa weak reference. Three-grid thermal Bose number matches independent p integral to relative <5×10⁻¹²; no rest-cell fake condensate.

## Critical negative and next action

At synthetic dt1e7, N8/16/32/64 net scalar changes ΔNφ=−7.59298/−7.28080/−6.85008/−6.42321×10⁻⁵. N32→64 still differs ~6.2%, N8→64 ~18.2%. No physical mode/time/p/c/UV convergence, and prior pmax2 missed ~65.6% of physical IR-curved source versus pmax16. Discontinuous P1 doubles endpoint degrees of freedom and may introduce discontinuity/threshold alias at energy interfaces. No general theorem guarantees positivity/implicit existence for every dt or proves physical BEC; zero-mode condensate requires separately defined δ³(p) dynamics. Full shared χ1χ2–SM/e±γν/Zprime+φ/FRW absent.

**NEXT_ACTION:** verify new real CI; compare 3+ parent pmax4/8/12/16 and parent angular12/24/48 while fixing physical mode support, compare continuous vs discontinuous energy hats at equal DOF, test nonthermal curved φ η and independent energy-interface traces, q/time/IR. Only then couple full χ1χ2 φ and shared finite SM plasma/FRW. No physical Γelastic, ΩDM/BEC, frozen fφ(q) or CLASS-ready assertion.
