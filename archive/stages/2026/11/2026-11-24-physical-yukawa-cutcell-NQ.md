# 2026-11-24 — physical Yukawa shell-crossing weak N/Q

Research stage: independent physically normalized positive weak daughter-bin integration for prescribed thermal χ1↔χ2+φ sources at masses 0.1/0.02/86e−9 GeV and Yukawa y=.001. No hidden bath, no physical nonlinear relic claim.

## Independent result

Cut each physical φ energy cell by exact boosted shell boundaries, split the *parent radial domain* at the exact shell/edge opening and closing roots, and evaluate actual Pauli/Bose occupation at each daughter event energy. The positive same-event measure yields signed Rφ, Qφ, Q1, Q2 and pointwise on-shell energy closure. A second independently derived log/dilog thermal daughter primitive followed by adaptive p quadrature agrees at 3 UV cutoffs, 3 parent/angle/φ resolutions and forward/inverse chemical tests.

At pmax=2 GeV, φ N16 first-cell Rphi=−2.762149360447387e−12 GeV4, Qphi=−1.227799761951597e−14 GeV5; correct physical total Rphi=−1.295453955491667e−11. The previously falsely zeroed first cell is resolved. At medium resolution nrad10 per parent segment ×16 log angles ×16 scalar cells, all four weak moments have relative errors ≲3.1e−12 against independent analytical 1D benchmark across P2/8/16 and three non-equilibrium cases. Physical common-T chemical equilibrium gives source near floating roundoff. All five focused local tests passed in 4.79 seconds.

**Negative result:** depositing the first-cell transferred number at its *physical DOS* energy mean W/V overstates its collision-weighted scalar energy-transfer magnitude by **45.16%** in the P2,N16 thermal-inverse state. Correct physical DOS capacities do not license point-mean strong moment closure. An unsplit 48-node parent quadrature still has ~2.92e−5 first-cell relative bias even if its total rate appears converged.

## Provenance and reproduction

- [Scientific derivation and numerical table](../../../../kinetics/PHI-BOSE-PHYSICAL-CUTCELL-NQ-AUDIT-2026-11-24.md)
- [Standalone source](../../../../kinetics/phi_bose_physical_cutcell_weak.py)
- [Focused tests](../../../../tests/test_phi_bose_physical_cutcell_weak.py)
- [Machine-readable values](../../../../kinetics/results/phi-bose-physical-cutcell-nq-2026-11-24.json)
- Execute: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_physical_cutcell_weak.py
- Execute: OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_physical_cutcell_weak

**Scope:** frozen thermal weak source, not a physical nonlinear strong solver, independent state capacity graph, 86-eV condensate, FRW relic Ω, CLASS or galaxy likelihood. New public Scientific CI result must be checked. Keep PR #1 draft and main unchanged.

**NEXT_ACTION:** construct a finite-capacity reversible on-shell strong operator using a realizable energy-resolving φ basis and common χ1/χ2/φ event phase measure; require exact N and all three Q weak limits, no spurious graph charges, finite-step Pauli/Bose bounds and entropy. If impossible under a specified discretization, retain the negative result and test alternative convex moment closure.
