# Constructive exact-energy Bose/Pauli resonance graph: three invariants, positive implicit microevents, but not continuous physics (2026-11-23)

## Motivation and precise scope

The [one-η energy/equilibrium no-go](PHI-BOSE-ONE-DOF-NO-GO-AUDIT.md) rejects cell-constant inverse fugacity as an exact continuum energy+thermal closure; the [earlier DG-P1](PHI-BOSE-DG-TANGENT-COUNTEREXAMPLE-AUDIT.md) rejects unqualified positivity of two-hat nonlinear moment inversion near ηφ=0. An independently testable intermediate construction uses **exact-resonant event energy states**. Its purpose is to prove that the desired algebraic requirements are simultaneously achievable on a *deliberately restricted physical on-shell graph*, then identify missing continuum/infrared convergence.

Use mχ1=.1, mχ2=.02, mφ=86×10⁻⁹ GeV and a synthetic energy lattice (GeV):
E1(i)=.100+.010i, E2(j)=.052+.010j, Eφ(k)=.048+.010k, with all reactions i=j+k, 0≤k≤i≤N. Each collision conserves E1=E2+Eφ **exactly**. For every event the boosted two-body energy Eφ=a−b cos θ* gives an allowed |cos θ*|≤1, verified explicitly; the at-rest first parent uses the exact shell within 10⁻¹⁰ GeV. A graph with N=2,4,8,16 has 6/15/45/153 physical shell events.

## Stoichiometric connectivity — no extra frozen linear charges on these tested graphs

Let the event incidence column be (−eχ1,i,+eχ2,j,+eφ,k). Its three left null vectors are the two number charges N1+N2, N1+Nφ and true discrete shell energy ∑Es Ns. On tested N=2,4,8,16, numerical ranks are 6,12,24,48 for 9,15,27,51 modes respectively, i.e. **left nullity exactly 3**. Max annihilation residual of the explicit invariants is ≤5.6×10⁻¹⁷. This matters because older adaptive collision graphs conserved energy and entropy while hiding false extra invariants. It is a *finite event network* statement, not a claim about the continuum collision operator's spectral gap or physical ergodicity.

## Quantum fixed point, admissible microstep, and time evolution

At common T=.7 and µ1=µ2+µφ with µ2=µφ=−.01 GeV, **every event** satisfies the exact FD/BE identity F=R; relative residual ≤6.7×10⁻¹⁶ at N≤16. For positive reaction extent ξ and positive mode capacities w, updated occupations are f1−ξ/w1, f2+ξ/w2, b+ξ/wφ. The intersection of Pauli/Bose constraints gives a *closed physical bracket* for ξ. The implicit equation ξ=dt λ(Fnew−Rnew) has at least one root by opposite signs at the two bracket ends: at the upper end F=0, at the lower end R=0. Bisection provides a positive occupation step **without clipping**. For any positive capacities, stoichiometric number and shell energy are preserved; physical discrete entropy S=∑w s(f,b) obeys ΔS≥dt λ(Fnew−Rnew)log(Fnew/Rnew)≥0 by concavity, provided the root is interior.

Single event f1=.5,f2=.2,b=.7, dt3, synthetic λ1, unit capacities: ξ=.207429193749, final (.292571,.407429,.907429), ΔS=.2549386 ≥ lower .04865635 and residual ~6×10⁻¹⁷. This is an algebraic reference only.

A connected N8 network with 45 microevents and synthetic λ=.02 was evolved by consecutive bounded implicit microsteps. Across 1/2/4/8 subdivisions of total illustrative time 2, Nφ changes by +.0491024 / +.0470730 / +.0463833 / +.0461727, with positive entropy and full global charge/energy errors ≤2.85×10⁻¹⁴. **Conservation does not imply time convergence:** the 1-vs-8-step difference is ≈6.35%. The tests cover two charges, energy, thermal fixed point, positivity, H and finite-step time sensitivity.

## Critical missing physics (do not promote to a cosmological solution)

This is a successful **discrete mathematical existence example**, not the actual physical φ production operator. Unit phase-space capacities and λ=.02 are arbitrary; the physical p²dp, spin/statistical degeneracies, Lorentz invariant decay weights, and continuum angular/angular shell integrals have **not** been assigned. In particular the lowest φ lattice energy is .048 GeV, whereas the true rest mass is 8.6×10⁻⁸ GeV: the ultralow-Eφ physical Bose tail and any q→0 Bose/δ³(p) mode are completely missing. Even a correct left nullspace and microreversible H theorem do not establish weak-rate convergence or late freezeout. Comparing to [exact mapped physical source](PHI-BOSE-DAUGHTER-LOG-ANGLE-AUDIT.md) requires physical capacities and events spanning true IR/UV, plus 3+ grids and independent bounds.

Reproduce with `OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_exact_shell_graph.py` or `python -m kinetics.phi_bose_exact_shell_graph`; four focused tests passed locally (~0.11 s). [Source](phi_bose_exact_shell_graph.py), [compact JSON](results/phi-bose-exact-shell-graph-2026-11-23.json). New GitHub CI after commit must be confirmed separately.

**NEXT_ACTION:** introduce physical positive on-shell phase-space capacities and event rates with low-Eφ shells (possibly adaptive shared resonant lattice), compare the reconstructed φ source and energy to the exact boosted-log-angle weak integral on 3+ radial/angle/IR grids; audit left nullspace after physical truncation and actual finite-step convergence. Only then consider full SM/FRW/86-eV observational claims.