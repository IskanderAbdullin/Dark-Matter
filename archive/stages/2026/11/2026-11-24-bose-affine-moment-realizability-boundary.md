# 2026-11-24 — exact Bose affine moment-face exclusion and positive non-affine counterexample

**Research result:** With independently converged physical χ1↔χ2+φ Yukawa weak first-cell Rphi/Qphi, the frozen moment ray M(t)=M(0)+t C crosses the **physical DOS integrated** DG-P1 *affine inverse fugacity* admissible face ηL=0 at P=2 GeV, scalar N=24: t*=174.23628158 GeV^-1; N32: t*=26.38466116 GeV^-1; N40: t*=149.52832484 GeV^-1. This improves prior **linear-tangent-only** diagnostics (N32 26.34743781) to exact finite moments.

At a fixed M_R, strict positive-definite Bose susceptibility implies ∂M_L/∂η_L|M_R=−(G_LL−G_LR²/G_RR)<0. Hence η_L=0 maximizes M_L in the affine η family; a strictly positive gap M_L(target)−M_L(max boundary) is **proof of exclusion** of an affine positive inverse-fugacity state. At 1.1t*, gaps are +4.16867e−13 / +2.34633e−13 / +1.50277e−13 GeV³ (N24/32/40), while for 0.8t* all three gaps are negative. N32 critical endpoint ηR*=0.0021512311069491253. Independent 50-digit mpmath quadrature confirms t*=26.3846611641.

**Counterexample for physical interpretation:** b_new=b_T(E)*(1+aLhL+aRhR) with positive physical DOS moments and coefficients from a 2×2 positive Gram system matches both *excluded* moments exactly and stays strictly positive: N32 at 2t* has aL=+1.51452e−4,aR=−3.29408e−4, i.e. both endpoint factors >0.99967. All N24/32/40 positive even at 10t*. Thus the barrier is specific to the numerical affine inverse-fugacity manifold, NOT nonnegative physical phase-space, and is not evidence of BEC, physical time blowup or a solved full Yukawa strong kernel.

At P=8/16 GeV and N16/24/32/40 local ηL tangent is positive; P=2,N48/64 first cell is truly forbidden by hard parent cutoff. No universal grid-independent IR physical claim.

[Derivation](../../../../kinetics/PHI-BOSE-MOMENT-REALISABILITY-BOUNDARY-AUDIT-2026-11-24.md), [source](../../../../kinetics/phi_bose_moment_realisability_boundary.py), [five tests](../../../../tests/test_phi_bose_moment_realisability_boundary.py), [data](../../../../kinetics/results/phi-bose-moment-realisability-boundary-2026-11-24.json). Five local focused pytest checks passed in 7.85s; public CI to verify independently.

**NEXT_ACTION:** directly evolve the true non-affine scalar modes in a *conditional fixed-fermion bath* with exact positive Yukawa on-shell scalar DOS, then compare finite-time moments to excluded affine-eta domain; only then advance to self-consistent shared-χ1/χ2/φ physical strong collision solver and FRW.
