# Stage 2026-11-22 — moving on-shell χ1↔χ2+real φ Bose/Pauli weak source

- **ID:** DM-KINETICS-20261122-MOVING-YUKAWA-BOSE-PAULI
- **Active session UTC start:** 2026-11-22T14:52:42Z. Prototype four tests rerun successfully at 15:34:11Z. Precise per-stage start/finish not separately measured.
- **Preceding branch:** [1424760e9a811b15288d7993a4791fff7d10d9b8](https://github.com/IskanderAbdullin/Dark-Matter/commit/1424760e9a811b15288d7993a4791fff7d10d9b8), public research/kinetic-equilibrium-audit-20261121 draft PR #1; new stage CI to verify.
- **Scope/router:** [GOAL](../../../../GOAL.md), [PLAN](../../../../PLAN.md), [kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md), [model](../../../../kinetics/MICROPHYSICS-MODEL.md), [gauge constraint](../../../../kinetics/REAL-PHI-YUKAWA-GAUGE-CHARGE-AUDIT.md), [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Scientific question and independent tests

The shared-lepton χ1/χ2 strong 2→2 solver still initializes χ2 synthetically, not via real φ Bose Yukawa production. Advance actual microscopic channel by deriving **moving-mother Lorentz-boosted** χ1↔χ2+φ number/energy Pauli+Bose WEAK operator, independently testing the existing scalar Yukawa width, threshold, full daughter energies, nonzero µ chemical equilibrium, both directions and physically meaningful sign reversal of net φ. It is deliberately **not yet** a coupled three-species strong/time cosmic rate.

## Physics and method

Real neutral φ Yukawa L_Y=−y φ(barχ2 χ1+h.c.), scalar trace |M|²avg=y²[(m1+m2)²−mφ²] and Γrest=|M|²avg sqrt λ/(16πm1³) for m1>m2+mφ. On-shell parent with p has E1=sqrt(m1²+p²), daughters boosted E2=(E1/m1)E2star+(p/m1)pstar c, Eφ=E1−E2, isotropic rest decay dc/2 and dilated Γ=Γrest m1/E1. Exact quantum J=f1(1−f2)(1+bφ)−(1−f1)f2 bφ. Source stoichiometry dotN1=−dotN2=−dotNφ, Q1+Q2+Qφ=0, entropy J ln(F/R)≥0. For bφ=0 inverse exactly zero; formal Sdot at b=0 is +∞ and code does not clip it. Chemical fixed point µ1=µ2+µφ; µφ need not be zero in the isolated closed channel.

[Script](../../../../kinetics/moving_yukawa_bose_pauli_weak.py), [4 focused tests](../../../../tests/test_moving_yukawa_bose_pauli_weak.py), [report](../../../../kinetics/MOVING-YUKAWA-BOSE-PAULI-WEAK-AUDIT.md), [JSON evidence](../../../../kinetics/results/moving-yukawa-bose-pauli-weak-2026-11-22.json).
Reproduction: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_moving_yukawa_bose_pauli_weak.py; python -m kinetics.moving_yukawa_bose_pauli_weak.

## Reproducible numeric result (synthetic, not 86eV mass)

Test m1=.1,m2=.02,mφ=.001,y=.001, Tφ=T2=1, parent FD T1=.9,1,1.1, pmax12. γrest=2.74969609751892e−9. With Gauss radial n8/16/32/48 and angular n12, net Rφ(n48) for T1=.9 is −2.65824136268e−11, T1=1 ≈1.71e−26, T1=1.1 +2.73810663852e−11. Thus reverse reaction actually destroys φ for χ1 colder than daughters while H remains positive (T1=.9: 7.419917e−12; T1=1.1: 6.934143e−12). Independent nonzero µ thermal condition at T=.8, µ1=.03,µ2=.05,µφ=−.02 yields Rφ≈4.88e−27, entropy≈1.21e−41; event four-energy defect ≤2.3e−15 and summed heat defect ~1e−26. Exact massless-daughter Γ=y²m1/(16π), and threshold m1≤m2+mφ closes reaction. Four local tests passed after fixing a vectorized array-angle validation error (initial tests failed for an array boolean comparison; fixed the actual input validation, not the physics tests).

## Limitation and next step

All occupations are **prescribed** for this WEAK integral, and real φ is not dynamically evolved/zero-mode treated; this cannot be combined with the old 5D shared-lepton strong time until a consistent finite-grid three-species Pauli/Bose **strong** event graph is derived with exact physical energy/number and entropy. mφ=.001 is **not** the 86eV benchmark, g1=4 numeric convention, parent p Gauss at these smooth thermal shapes does not prove global nonthermal/IR daughter spectral convergence. No actual SM+mediator finite bath or FRW, ΩDM/Γelastic/BEC/CLASS-ready. **Next P0:** energy-matched fully positive χ1/χ2/φ reaction graph with accurate phase-space measures, Fermi/Bose bounds, full N/E/H and independent q/dt/IR/cutoff; integrate with shared-lepton strong source only after those gates. Archive and checkpoint remain permanent. User's 55-minute run target does not imply background work.
