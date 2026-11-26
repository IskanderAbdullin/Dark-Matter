# 2026-11-24 — converged Bose first-cell DG-P1 inverse-fugacity tangent negative

After the independently converged physical parent-shell cut-cell 3-species weak N/Q stage, the source of the scalar first cell was projected onto the two local entropy-dual affine hats using exact physical Rphi and Qphi. A separate physical p²dp Bose susceptibility gave unique dη/dt=−G^-1 C. At P=2 GeV, scalar N16 has dηleft=+5.70429065e−9 GeV, N24 −4.93472804e−10, N32 **−3.26407450e−9**, N40 −5.76416044e−10, and N48/64 exactly zero because no kinematically reachable first-cell daughter energy.

At N32 both physical scalar hat sources are negative, but the left inverse fugacity tangent is negative. The local linear ηrest positivity boundary estimate is 26.3474 GeV^-1, **not** an implicit-step no-existence or BEC theorem. This is stable over independent 6×8,10×16,16×32,32×48 parent/angle quadrature and the susceptibility agrees with its massless high-T asymptote to <1.5e-4 relative. Earlier angular-only N32 coefficient was ~7.34% higher magnitude, establishing the need for exact shell-edge splits.

Source [new script](../../../../kinetics/phi_bose_exact_cutcell_tangent.py), [four focused tests](../../../../tests/test_phi_bose_exact_cutcell_tangent.py), [primary data](../../../../kinetics/results/phi-bose-exact-cutcell-tangent-2026-11-24.json), [report](../../../../kinetics/PHI-BOSE-EXACT-CUTCELL-TANGENT-GATE-2026-11-24.md). Four tests passed locally in ~6.82s; public CI pending.

NEXT_ACTION: test constrained physical realizability in finite-time energy-resolving χ1/χ2/φ collision moments; prove Pauli/Bose admissibility, N/E/entropy, equilibrium and true weak limit. No artificial clipping, condensate or physical Ω claim.
