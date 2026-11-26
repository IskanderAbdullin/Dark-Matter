# Exact daughter-energy cut-cell moments: analytic weak-measure control (2026-11-23)

**Status: analytic kinematic identity and local independent numerical spot checks; NOT a physical strong finite-volume collision operator.** Parent on-shell boost with E_phi(c)=A+B c, c in [-1,1], A=gamma E_phi*, B=gamma beta p*, B>0. For daughter-energy bin [L,U), intersect [L,U) with [A-B,A+B]: l=max(L,A-B), h=min(U,A+B). If h<=l, both moments vanish. Otherwise

- angular probability P_bin=(h-l)/(2B);
- angular first energy moment M_phi,bin=(h²-l²)/(4B) = P_bin (h+l)/2;
- accompanying daughter-2 energy moment M_2,bin=E_1 P_bin - M_phi,bin (since E_2=E_1-E_phi pointwise);
- parent energy moment M_1,bin=E_1 P_bin.

For B=0 use the continuous limit: bin indicator of A, P=1 or 0, M_phi=A P. Across a complete partition covering [A-B,A+B], sum P=1, sum M_phi=A, sum M_2=E_1-A. Thus the **weak** number and three energy moments obey exact event-level kinematic energy accounting, including partially intersected daughter bins. These identities hold for angular-independent amplitude/occupation factor; if occupations depend on daughter energies, integrate the factor inside the bin, do not replace it with a midpoint value and claim exactness. Spin multiplicity, invariant parent measure, matrix-element normalization and shared physical daughter capacities are NOT supplied by these identities.

## Independent local check

Command: Python 3, SciPy integrate.quad over c in [-1,1] with split points at clipped bin edges; compare analytic P and M_phi to direct integration of 0.5*1_bin and 0.5*(A+B*c)*1_bin.

| A | B | [L,U) | analytic P | analytic Mphi | quadrature P | quadrature Mphi |
|---:|---:|---|---:|---:|---:|---:|
| 2 | 0.6 | [1.6,2.1) | 0.4166666666666667 | 0.7708333333333333 | 0.4166666666666667 | 0.7708333333333334 |
| 2 | 0.6 | [2.8,3.1) | 0 | 0 | 0 | 0 |
| 2 | 0.6 | [1,3) | 1.0000000000000002 | 2.0000000000000004 | 1 | 2 |
| 2 | 0.00001 | [1.999999,2.000002) | 0.14999999998766 | 0.30000007498198755 | 0.14999999998765998 | 0.30000007497531994 |

**Negative/limitation:** A positive angular cut-cell weak integral does not by itself furnish a reversible shared-event strong operator: daughter bin-averaged occupations, exact discrete energy, capacities, and H theorem cannot all be assumed. Numerical cancellation in (h²-l²) for very small B can be avoided by M_phi=P*(h+l)/2. For generic energy-dependent Pauli/Bose factors the above unweighted moments are only geometric components, not physical integrated rates.

**NEXT_ACTION:** implement independently normalized parent-radial × angular cut-cell integration with actual Yukawa gain/loss occupations, compare all three Q and N to continuous reference at >=3 parent/angle/daughter resolutions; only then derive shared positive phase-space capacities and strong nonlinear operator. No BEC/relic claim.
