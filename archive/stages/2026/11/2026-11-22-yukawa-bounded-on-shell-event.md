# Stage 2026-11-22 — bounded exact-energy χ1↔χ2+real φ Bose event (single mode triplet)

- **ID:** DM-KINETICS-20261122-YUKAWA-BOUNDED-ON-SHELL-EVENT
- **UTC:** full active session started 2026-11-22T14:52:42Z; five tests completed around 15:39:50Z; no claimed exact stage duration.
- **Starting branch:** [960980e5f59ee7878000a4b107e0cfad6f918c11](https://github.com/IskanderAbdullin/Dark-Matter/commit/960980e5f59ee7878000a4b107e0cfad6f918c11), public research/kinetic-equilibrium-audit-20261121 draft PR #1, stage CI pending.
- **Scientific routing:** [GOAL](../../../../GOAL.md), [PLAN](../../../../PLAN.md), [kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md), [Yukawa moving weak](../../../../kinetics/MOVING-YUKAWA-BOSE-PAULI-WEAK-AUDIT.md), [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Goal, test plan and on-shell strong micro-event

Before any physical finite-grid φ Bose time, test one EXACT energy-matched moving decay mode χ1↔χ2+real φ. Derive bounded implicit reaction extent x, requiring **no Fermi/Bose clipping**, exact Nχ1+Nχ2 and Nχ1+Nφ and total E, nondecreasing quantum entropy, distinct forward/inverse chemical equilibrium and rigorously closed threshold from parent. Test 250 independent extreme nonthermal cases. Alternative: unsafe cell interpolation and false entropy despite conserved energy (previous χℓ negative).

From boosted χ1 mother M1=.1,p=.3,c*=.2 to χ2 M2=.02 and φ Mφ=.001, E1=.31622776601683794, E2=.1932193769133251, Eφ=.12300838910351285. Choose arbitrary **positive test modal capacities** w1=.1,w2=.15,wφ=.2 and event W=.05. f1=a−x/w1,f2=d+x/w2,b=b0+x/wφ, quantum J=f1(1−f2)(1+b)−(1−f1)f2 b, implicit x=dt W J(x). Bracket x_lower=max[(a−1)w1,−d w2,−b0wφ], x_upper=min[a w1,(1−d)w2]; positive/negative x represent physical decay/inverse. Solve using Brent; finite entropy with scipy xlogy at empty b=0. This is NOT exact continuous-time integration, not real physical f(q) cell capacities, and not a global multi-event solver.

## Reproduction and results

[Code](../../../../kinetics/yukawa_on_shell_bounded_event.py), [5 tests](../../../../tests/test_yukawa_on_shell_bounded_event.py), [full audit](../../../../kinetics/YUKAWA-ON-SHELL-BOUNDED-EVENT-AUDIT.md), [machine numbers](../../../../kinetics/results/yukawa-on-shell-bounded-event-2026-11-22.json). Run OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_yukawa_on_shell_bounded_event.py. Five tests passed locally (~0.57s).

At dt10: initial (f1=.6,f2=.3,b=.1) gives decay x=+.0328142120, ΔS+.0659934336; inverse initial (.05,.8,2) x=−.0365123342, ΔS+.0605746330; empty Bose initial (.4,.1,0) grows b=.1423904 after x=+.0284780807, ΔS+.0958931627. Exact charge/energy defects ~1e−17. Independent nonzero chemical μ with μ1=μ2+μφ fixed point within1e−13. Deterministic seed20261122 stress across 250 Fermi/Bose/capacity/stiffness combinations produced no physical bound/H/charge failures.

**Negative diagnostic:** initial test used fixed abs conserved tolerance3e−12 but found rounded 3.64e−12 at quantity18858 (≈2e−16 relative). Corrected to explicit absolute5e−13 plus relative5e−14; does not modify model or hide a physical error. All three levels use **arbitrary w1/w2/wφ**, so do not claim continuum physical phase-space, realistic Bose IR or full χ2 production history. Only on-shell elementary microphysics verified.

## Next P0

Verify real GH Scientific CI, then construct **physical** finite-volume χ1,χ2,φ mode graph with energy-conserving exact kinematic reactions and actual cell capacities/interpolation proved entropy-compatible; compare 3+ q/angular/IR/cutoff and independent boosted moving weak, integrate into TWO χ+ONE shared lepton strong Pauli bath with real Yukawa replenishment and complete Zprime/SM/FRW. Physical mφ86eV, γelastic/H, ΩDM, BEC, frozen f(q), CLASS-ready remain unproven.
