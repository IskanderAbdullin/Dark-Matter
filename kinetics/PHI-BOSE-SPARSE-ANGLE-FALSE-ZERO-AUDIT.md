# Negative result: sparse Gauss-angle Yukawa events create a **false empty first Bose energy cell** (2026-11-23)

## Precise falsified inference

The earlier [DG-P1 cell entropy audit](PHI-BOSE-CELL-ENTROPY-AUDIT.md) correctly reports that, at scalar quadratic N=16, parent pmax=2 GeV, **its finite 12×12 Gauss event graph** gives zero source and zero inverse-fugacity tangent in the first scalar cell. This is a valid algebraic statement about that event graph, but it must **not** be interpreted as a physical on-shell kinematic absence. We now compare the actual continuous boosted decay support to the minimum energy of the sampled Gauss events **using the exact original fixture's grid definition**.

With \(m_1=.1\), \(m_2=.02\), \(m_\phi=86\) eV, quadratic scalar energy edges
\[
E_j=m_\phi+\left[\sqrt{P^2+m_1^2}+0.2-m_\phi\right](j/N)^2,
\]
the first cell is \([m_\phi,E_1]\). The continuous on-shell minimum over \(0\le p\le P\) follows the [exact rapidity inversion](PHI-BOSE-IR-ACCESSIBILITY-AUDIT.md). The original source samples only \(p_i=P(1+x_i)/2\) and rest-frame angles \(c_j\), where \(x_i,c_j\) are **12-point Legendre-Gauss roots**. Since no quadrature node reaches \(p=P\) or \(|c|=1\), the sampled minimum can be much larger than the true minimum.

## Direct independent numerical result

| Parent P (GeV) | True scalar minimum (GeV) | Minimum among 12×12 sampled events (GeV) | Sampled / true |
|---:|---:|---:|---:|
| 2 | 0.00119925094 | 0.00936825671 | 7.81 |
| 8 | 0.000299988288 | 0.00954935121 | 31.83 |
| 16 | 0.000149998548 | 0.00990482446 | 66.03 |

At the **previously discussed N16, P=2** the first scalar cell top is **0.00860359519 GeV**. It is **larger** than the true 0.00119925-GeV minimum but **smaller** than the sampled 0.00936826-GeV minimum. Therefore the first cell has **nonempty positive-measure physical on-shell support**, yet the 12×12 event graph deposits exactly zero there. This is a **quadrature false zero**, not a kinematic zero.

A genuine physical first-cell no-support example occurs for **N64, P=2**: first-cell top ≈0.000537805325 GeV, below the true minimum ≈0.00119925094 GeV. That cell is kinematically empty even with exact continuous angular integration. Other combinations can be classified using the supplied audit script rather than a generic rule.

The original DG-P1 algebraic conclusion survives unchanged: nonoverlapping susceptibility blocks prevent a first-cell inverse-fugacity tangent when **the discrete source is exactly zero**. But the stronger physical inference is **refuted**: zero discrete first-cell source does not imply zero continuous Yukawa collision source. Likewise the old N16 dt30 Newton acceptance is not a physical IR convergence certificate.

## Why this matters

- An event-graph nullspace can be caused by a hard parent cutoff, a truly forbidden kinematic region, **or insufficient quadrature resolution**. These are distinct mechanisms and must be separated.
- Increasing scalar cell resolution without increasing parent/angle resolution can make the first-cell source disappear spuriously. A green energy/charge/H theorem and implicit time solver cannot detect this omission.
- Near-critical Bose amplification and possible low-energy spectra are especially vulnerable to this numerical false zero. It cannot be repaired by inventing a condensate or clipping Bose inverse fugacity.

## Reproduction and decision

[Source](phi_bose_sparse_angle_false_zero.py), [focused tests](../tests/test_phi_bose_sparse_angle_false_zero.py). Run:

\`\`\`bash
OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_sparse_angle_false_zero
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_sparse_angle_false_zero.py
\`\`\`

The threshold calculations were independently reproduced in a local Python/NumPy 12-point Gauss script; public GitHub CI of the new tests must be verified separately. No full continuous first-cell integrated collision magnitude has yet been computed, so this is a **support counterexample**, not a quantified first-cell rate correction.

**Next decisive test:** replace fixed Gauss-angle sampling by **cut-cell adaptive angular integration** with exact on-shell daughter-energy interval intersections for each scalar cell; independently converge parent radial 12/24/48, angle 12/24/48 or exact per-cell k integration, scalar N=8/16/32/64, P=2/8/16. Compare each cell's physical rate and energy moment to independent 1D chemical source, then repeat nonzero-μ equilibrium, inverse-dominated and finite-time entropy tests. Physical strong finite-volume capacity, SM/FRW and 86-eV BEC/Ω/CLASS remain unverified.
