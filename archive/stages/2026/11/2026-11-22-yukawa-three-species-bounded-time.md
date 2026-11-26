# Stage 2026-11-22 — moving χ1↔χ2+real φ nonlinear Fermi/Bose finite-volume kinetics

- **ID:** DM-KINETICS-20261122-MOVING-YUKAWA-THREE-SPECIES-TIME
- **Active session START_TIME_UTC:** 2026-11-22T15:38:09Z, measured with container date. Per-stage ending UTC not separately measured; do not infer 55 minutes from action count.
- **Repository/branch:** public IskanderAbdullin/Dark-Matter, research/kinetic-equilibrium-audit-20261121, draft PR #1. **Prior HEAD:** [960980e5f59ee7878000a4b107e0cfad6f918c11](https://github.com/IskanderAbdullin/Dark-Matter/commit/960980e5f59ee7878000a4b107e0cfad6f918c11), prior [Scientific CI #45039365150](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45039365150) success on independent moving Yukawa WEAK.
- **Context:** [repository GOAL](../../../../GOAL.md), [PLAN](../../../../PLAN.md), [AGENTS](../../../../AGENTS.md), [kinetics specification](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md), [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md), [current state](../../../../kinetics/research-control/STATE.json), [next action](../../../../kinetics/research-control/NEXT_ACTION.json).

## A priori question and physical restrictions

Prior moving χ1→χ2+real φ on-shell Yukawa collision was only a weak integral evaluated at prescribed Fermi χ1/χ2 and Bose φ functions. Test actual **three independent evolving occupation spectra** under a common positive Lorentz-boosted moving-parent microscopic event with correct Bose stimulation, daughter Pauli and inverse χ2+φ→χ1. Demand N1+N2 and N1+Nφ, exact exchanged energy and nonnegative actual entropy at EVERY accepted time step and Fermi 0<f<1, Bose b>0, not merely global weak Q. Independent thermal weak/strong match for 3 energy grids and nonzero separate μ chemical equilibrium; preserve negative outcomes at coarse grid and in reverse direction.

## Exact on-shell event and discrete H proof

Use original [physical scalar width and boosted daughters](../../../../kinetics/moving_yukawa_bose_pauli_weak.py), m1=.1,m2=.02,mφ=.001,y=.001 in **arbitrary energy units**, not 86eV physical relic. Per physical moving mother momentum p, E1=√(p²+m1²), E2+Eφ=E1, and isotropic rest-frame angle c, event measure g1 p²dp/(2π²) Γrest m1/E1 dc/2. Quantum forward F=f1(1−f2)(1+bφ), inverse R=(1−f1)f2 bφ, J=F−R.

For each species energy grid, barycentric α(E) exactly reproduces unity and E, so paired −J,+J,+J deposition preserves both stoichiometric null combinations and **total physical energy**. The same interpolation is used in the event occupations: fermion logit h_i=ln[f_i/(1−f_i)], Bose η_i=ln[(1+b_i)/b_i]. At each event h(E)=Σα_i(E)h_i and η(E)=Σα_i(E)η_i; Fermi f=expit(h), Bose b=1/expm1(η). Then the physical nodal S with radial density-of-states w=g/(6π²)[(U²−m²)^(3/2)−(L²−m²)^(3/2)] has

    dS/dt = -Σ_i dotN1_i h1_i -Σ_j dotN2_j h2_j +Σ_k dotNphi_k ηphi_k
           = ∫dW (F−R)ln(F/R) >=0.

This is a **semi-discrete source H theorem**, not a rigorous finite-step global error guarantee. Finite time uses adaptive full-step/two half-Euler spectral difference (f1,f2,log1p b) plus hard reject on any Fermi/Bose bound or actual entropy decrease, no clipping or thermal reset. First order. Zero-Bose endpoint b=0 intentionally **fails closed** pending a proper physical limiting method; no BEC or zero-mode physics.

## Real test evidence and failures

- [Full runnable three-species source](../../../../kinetics/yukawa_bose_pauli_bounded_time.py), [8 focused pytest functions / 12 tests with prior weak suite](../../../../tests/test_yukawa_bose_pauli_bounded_time.py), [scientific audit](../../../../kinetics/YUKAWA-THREE-SPECIES-BOUNDED-TIME-AUDIT.md).
- [Six synthetic physical-cell evolution runs](../../../../kinetics/results/yukawa-three-species-pauli-bose-time-2026-11-22.json) and [three local tolerance runs](../../../../kinetics/results/yukawa-bose-pauli-temporal-tolerance-2026-11-22.json).
- Exact commands: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_yukawa_bose_pauli_bounded_time.py tests/test_moving_yukawa_bose_pauli_weak.py → **12 passed locally, 2.21s**. OPENBLAS_NUM_THREADS=1 python -m kinetics.yukawa_bose_pauli_bounded_time.
- Deterministic Gauss mother p nr12, rest angle12, pmax12; grids16/32/64 energies each for χ1,χ2,φ; g1=g2=4/gφ=1 **numeric historical convention**, not validated physical relativistic SM/gauge counting. y=.001, toy end time5e7, initial dt2e7, atol1e−6.

At identical thermal FD/BE f( E) (T1=1.1, T2=Tφ=1), full weak R,Qχ1,Qχ2,Qφ,S equals independently implemented [prior moving-Yukawa weak](../../../../kinetics/moving_yukawa_bose_pauli_weak.py) to relative 1e−10 for each bins16/32/64. T=.8, μ1=.03, μ2=.05, μφ=−.02 exact separate chemical fixed point. In six actual nonlinear histories (T1=.9 or1.1), all both conserved charges and total E residual ≤2e−15, all accepted entropy increments >=0, f1,f2∈(0,1), bφ>0. At bins32,T1=1.1: 127 accepted, 48 rejected, ΔNχ2=ΔNφ≈+0.00126852, actual entropy +0.000324541; T1=.9 gives NET φ destruction ΔNφ≈−0.00119416 and entropy +0.000337165. These are synthetic toy kinetic results, not cosmic evidence.

**Genuine NEGATIVE:** bins16→64 moves |ΔNφ| by ~8.9% for T1=.9 and ~2.6% for T1=1.1, despite essentially exact conservation and H theorem. At fixed bins32 T1=1.1, reducing local atol 1e−5→1e−6→1e−7 increases accepted/rejected 40/15→127/48→396/140, while ΔNφ shifts +.0012691535→+.0012685176→+.0012683177 (first-order bias). No convergence certificate for Bose low-energy IR or radial/angular/UV and no physical 86eV benchmark.

## Conclusion, peer-scientific context, continuation

This establishes a *scoped* coupled evolving **χ1,χ2,real φ** finite-volume Pauli/Bose decay pilot with number/energy/entropy and actual occupations, now beyond the prior prescribed-source weak test. See [Yin 2023](https://arxiv.org/abs/2301.08735) and [Sakurai & Yin 2025](https://arxiv.org/abs/2410.18968) for the Bose-enhancement motivation, not borrowed code. It does not have the one-shared-lepton χ1/χ2 vector elastic system connected; no SM e±γν+Z′ finite bath, FRW H/g*, IR zero-mode, physical masses/thermal constraints, late freezeout, ΩDM or CLASS-ready. No scientific discovery claimed.

**NEXT P0:** couple one shared lepton elastic generator with this exact 3-species Yukawa generator into a 4-population χ1/χ2/φ/ℓ semi-discrete positive-entropy reaction network with correct nullspace (total dark fermion charge, conserved N1+Nφ under pure channels, separate lepton number) and energy, step-by-step Pauli/Bose. Then investigate ≥3 bins/radial/angle/time/cutoff/IR and actual physical parameters. No merge main/force; permanent negative data retained.
