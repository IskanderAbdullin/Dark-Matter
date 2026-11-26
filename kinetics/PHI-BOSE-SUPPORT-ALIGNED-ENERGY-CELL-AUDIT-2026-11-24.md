# Exact scalar support-edge-aligned Bose DG-P1: eliminate a false negative IR tangent (2026-11-24 UTC)

## Why this test matters

The independently converged physical Yukawa weak cut-cell source at parent P=2 GeV gave a **negative left inverse-fugacity tangent** in the original scalar first-bin DG-P1 (N=24/32/40), despite both scalar hat collision moments being negative. A physically positive exact scalar-mode finite-time solution with fixed FD1/FD2 was shown to leave the **positive affine η first-cell moment family** at N32 t≈29.023 GeV^-1. This is a numerical representation obstruction, not BEC.

**New geometric diagnosis:** the original quadratic first energy bin [mφ,U_N] actually *straddles* the exact physical on-shell scalar support floor

    k_floor(P2) = 0.0011992509375794726 GeV,

while mφ=86e−9 GeV. Below k_floor, the isolated finite-P=2 on-shell 1↔2 Yukawa collision kernel is identically zero, regardless of occupations. Above it, inverse reactions deplete the scalar spectrum. A single affine entropy-dual η(E) across this internal **physical support kink** attempts to represent two different dynamics simultaneously and may falsely require η(mφ)<0. This diagnosis is falsifiable: split the first cell at k_floor and recompute the physical two-moment susceptibility and finite-time exact-mode moments.

## Method: split exactly at an ON-SHELL boundary, not a made-up condensate cell

For each original N=24/32/40 first scalar cell [mφ,U], replace it with two **nonoverlapping** cells [mφ,k_floor] and [k_floor,U]. The lower cell receives **exactly zero** on-shell collision events and retains its initial thermal Bose occupation in the *fixed-fermion, finite-parent-only* scalar modal equation. In the upper supported cell, the physical weak number and scalar heat Rφ,Qφ are exactly the original first-bin sources (no events below k_floor).

For any energy segment [L,U], define affine local hats hL=(U−E)/(U−L), hR=(E−L)/(U−L); source projections C=( (U R−Q)/(U−L), (Q−L R)/(U−L)). Independently integrate its positive physical massive Bose susceptibility G_ab=(1/2π²)∫p²dp h_a h_b b(1+b), and take dη/dt=−G^-1 C. Crucially G uses the **supported momentum subinterval** p∈[sqrt(L²−mφ²),sqrt(U²−mφ²)] rather than the old full first-bin DOS.

## Decisive tangent sign reversal with **unchanged physical weak R/Q**

The original physical parent/angle shell-crossover weak fluxes are unchanged; only the spectral basis is aligned with a real support boundary.

| Scalar first-bin N at P2 | Original dη_left/dt (GeV) | Aligned active-cell dη_left/dt (GeV) | Aligned active-cell dη_right/dt (GeV) |
|---:|---:|---:|---:|
| 24 | **−4.93472804e−10** | **+7.45138036e−9** | +1.11498398e−8 |
| 32 | **−3.26407450e−9** | **+2.94300213e−9** | +1.26229221e−8 |
| 40 | **−5.76416044e−10** | **+2.07608769e−10** | +4.95693292e−9 |

This removes the previously diagnosed local Bose positivity tangent obstruction in **all three tested grids**, without changing the event rate, inventing IR emission, clipping inverse fugacity or adding a condensate. The inactive subcell has exactly zero source and its η remains unchanged under the restricted 1↔2 channel. All tested local G matrices are strictly positive definite and satisfy G·dη+C=0 to floating precision.

## Beyond tangent: genuine exact positive scalar finite-time TWO-MOMENT fits

Use the previously derived **exact scalar-mode semigroup** with physically normalized Yukawa S(k),D(k) and frozen FD fermions: b(k,t)=b0 exp(−Dt)+S*t*exprel(−Dt). Integrate its actual two physical p²dp moments in the supported subcell [k_floor,U] with positive daughter-k quadrature, not a frozen linear ray. Then solve for a positive affine inverse fugacity η(E)=ηL hL+ηR hR with ηL,ηR>0 that matches both physical moments. Positivity is enforced by ηL=L exp(xL), ηR=U exp(xR); matching two moments alone is **not** a claim of matching the full non-affine spectrum pointwise.

Compare actual fixed-bath time t to the earlier **original unsplit frozen-ray** critical time t* for each N:

| N | Earlier unsplit frozen t* (GeV^-1) | Supported-cell ηL,ηR at t=2t* (GeV) | Supported-cell ηL,ηR at t=100t* (GeV) |
|---:|---:|---|---|
| 24 | 174.23628158 | (0.001201850024,0.003827754925) | (0.001335453438,0.004022306364) |
| 32 | 26.38466116 | **(0.001199406236,0.002151629504)** | **(0.001207011382,0.002184526185)** |
| 40 | 149.52832484 | (0.001199312897,0.001378130573) | (0.001202035407,0.001452533712) |

For all three grids, positive fitted ηL and ηR exist at **1.1,2,10,100 times** the previous old frozen-ray critical time, with normalized two-moment residual ≤5e−12. At N32 t=2t*, the independently evaluated fitted eta values are 0.00119940623606 and 0.00215162950393 GeV, versus the allowed lower supported-edge etaL0≈0.00119925093758 and upper edge etaR0≈0.00215096329829 GeV. At 100t* the true active-cell moments converge across independent (scalar Gauss,parent Gauss) orders 12×12,24×16,48×24,96×48,128×64 to <2e−15 fractional difference for N32; this is **numerical evidence of reproducibility**, not a proof of all-time representability.

The original *unsplit* N32 positive affine inverse-fugacity family remains mathematically unable to represent the exact positive scalar modal moments even at t=1.1t* (positive fixed-MR gap previously proven), and that negative result is not erased. The two cells have different bases and capacities; the split resolves a representation failure rather than making the original negative evidence disappear.

## Scope and reproducibility

Source: [support-aligned finite-moment and tangent audit](phi_bose_support_aligned_moment_closure.py). Tests: [three targeted regressions](../tests/test_phi_bose_support_aligned_moment_closure.py). Commands:
- OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_support_aligned_moment_closure.py
- OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_support_aligned_moment_closure

Three local focused tests passed in 4.61s (physical support split, tangent sign reversal, finite-time physical two-moment fits across N24/32/40 and 1.1–100 times old boundary, N32 exact eta numerical regression). Public Scientific CI separately verify.

**Scientific boundary:** this is only the first scalar energy cell and fixed Fermi bath with hard parent cutoff P=2; *full pointwise spectral shape* is not matched by the two-moment fit, and other cells may still become inadmissible. If parent UV cutoff changes, physical support floor moves, and a static split cannot be assumed to solve global IR time evolution. Other channels, SM/FRW/redshift and a real zero mode may couple the presently inactive region. The work does **not** establish a fully coupled Pauli/Bose conservative strong finite-time method, BEC or actual relic Ω.

**NEXT_ACTION:** make physical on-shell support discontinuities first-class adaptive energy cell boundaries (including time-dependent q/FRW when implemented), and test the *entire* three-species positive finite-step energy-resolving solver with correct shared phase capacities, both charges, Q ledger, entropy, gap/nullspace and cutoff/time convergence. Do not use old unsplit DG-P1 false negative as physical evidence of condensate.
