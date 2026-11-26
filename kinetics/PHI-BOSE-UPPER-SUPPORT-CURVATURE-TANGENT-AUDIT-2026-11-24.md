# Physical upper-support curvature refutes global DG-P1 Bose tangent positivity, even after support alignment (2026-11-24 UTC)

## Why the earlier IR repair is not sufficient

The independently converged P2 Yukawa χ1↔χ2+φ scalar weak benchmark (m1=.1, m2=.02, mφ=86e−9 GeV, y=.001, FD1 at T1=.9, FD2/BEφ at T2=Tφ=1, zero μ) has a **true scalar on-shell support interval** [kmin(P),kmax(P)]. Splitting the *first* Bose DG-P1 energy cell at kmin repaired a severe negative inverse-fugacity tangent. It does **not** guarantee a nonnegative inverse-eta tangent at every endpoint of the entire spectrum.

At P=2 GeV, the scalar upper support is

    K=kmax(P)=1.9211992509359777 GeV.

The original Nφ=32 quadratic grid's last active scalar bin is j=29, [L=1.8088878938628088,K], with no physical collision events above K. The physical χ1↔χ2+φ reverse reaction is more probable than forward for these thermal state temperatures, so **every true physical on-shell scalar mode** is depleted and its pointwise inverse fugacity rate is nonnegative, vanishing exactly at K. Yet the positive physical DOS-weighted affine-η two-hat moment projection yields a **negative η derivative at the right endpoint**.

## Modewise sign and exact projection mechanism

Using \(f_1(E;T_1=.9)\), \(f_2(E;T_2=1)\), \(b_\phi(k;T_\phi=1)\), and μ=0, the exact quantum odds give

    F/R = exp[ E1(1/T2−1/T1) ] = exp[−E1/9] < 1

for every actual on-shell event (where F/R is defined). Thus J=F−R<0 at every event and the true scalar modal derivative is \(\dot b(k)<0\) everywhere in its open on-shell support. Therefore the true inverse-fugacity derivative

    g(k)=dη(k)/dt = −(db/dt)/[b(k)(1+b(k))]

is **strictly positive** inside physical support, while g(K)=0 because the inverse parent event interval shrinks to measure zero. This is a continuum physical sign property; no condensation is involved.

For an affine DG-P1 two-moment Bose entropy-dual basis on a physical scalar energy subcell [K−h,K], the local susceptibility is \(G_{ab}=\int_{K-h}^K dk\,D_\phi(k)b(1+b)h_a h_b\) and sources C_a=\(\int D_\phi h_a\dot b\,dk\). Solving \(\dot\eta=-G^{-1}C\) is **exactly** the weighted L² projection of the true positive g(k) onto affine functions. It need not preserve pointwise positivity of g at either endpoint even though the event rate and all physical moments are correct.

The independent scalar-mode Yukawa integral gives the local asymptotic expansion with δ=K−k:

    g(K−δ)=A δ + B δ² + O(δ³),
    A≈2.3823467443e−11, B≈7.8894549565e−12 >0

in natural units. This is obtained by a quadratic fit to g(δ)/δ at δ=.001–.006 GeV, with actual two-sided parent support and physical Bose phase DOS; it is independently checked against the moment projection. For a smooth positive physical susceptibility weight, the affine weighted L² projection of δ² onto span{1,δ} on [0,h] has at δ=0 the leading value **−h²/6** (since for uniform weight P[δ²]=hδ−h²/6; nonconstant smooth weight first changes higher orders). Hence the projected **right endpoint** inverse-eta derivative is

    [dotη_R]_DG-P1 = −(B/6) h² + O(h³) < 0

for sufficiently narrow supported cells, despite **true g(K)=0 and true g(k)>0**. This is an explicit mathematical explanation of the persistent numerical tangent sign defect.

## Reproducible numerical negative result

Integrate the scalar source Rphi and Qphi over only the physical supported [K−h,K] interval using independent *scalar energy-first, parent inverse-support second* quadrature, not the parent-first angle graph. Compute physical massive scalar DOS Bose susceptibility separately by p²dp. At Nφ32, its original last active supported width is h0=0.1123113570731689 GeV.

| Subdivide h0 by | Width h (GeV) | Computed \(\dotη_{R}\) (GeV) | Ratio to \(-Bh²/6\) |
|---:|---:|---:|---:|
| 2 | 0.0561556785 | −4.26575444e−15 | 1.02876 |
| 4 | 0.0280778393 | −1.05132775e−15 | 1.01418 |
| 8 | 0.0140389196 | −2.60982519e−16 | 1.00704 |
| 16 | 0.0070194598 | −6.50168549e−17 | 1.00351 |
| 32 | 0.00350972991 | −1.62257651e−17 | **1.00176** |
| 64 | 0.00175486495 | −4.05289435e−18 | **1.00088** |

All right-node derivatives remain negative as the last supported interval is refined; their magnitude decreases approximately by factor four for each halving h, in agreement with the independent curvature asymptote. The left-node projected derivative remains positive. The entire original last-bin Rφ and Qφ computed scalar-first agree with independent exact parent-shell cut-cell weak values within 1e−11 relative.

**Negative result:** splitting the grid at both exact kmin and kmax does **not** force a positive DG-P1 inverse-eta derivative at all local endpoints. The first-cell support split is a significant conditional IR improvement, but no all-cell positivity theorem follows. The high-edge residual is numerically small and vanishes quadratically with refinement; its **sign** remains an exact projection problem, not a sign of physical stimulated Bose creation or a true cutoff-zero-mode.

## Tests and next step

Source: [phi_bose_upper_support_curvature.py](phi_bose_upper_support_curvature.py).
Focused tests: [test_phi_bose_upper_support_curvature.py](../tests/test_phi_bose_upper_support_curvature.py).
Run commands:
- OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_upper_support_curvature.py
- OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_upper_support_curvature

Three focused local tests passed in 1.20s; microscopic pointwise sign and curvature A/B, independent full last-bin scalar-first R/Q, subcell h and h² right-node scaling. Public CI separately verify. The simple asymptotic coefficient B is a *numerical fit* to the physical modal source; its sign and high-order dependence are supported by independent projection convergence, not claimed as a rigorous symbolic derivative formula for arbitrary model parameters.

**NEXT_ACTION:** a truly robust physical kinetic discretization must be more than DG-P1 support-edge alignment. Explore a positive shape-preserving/non-affine entropy-variable or direct occupancy representation with exact physical DOS, common on-shell microevent extents and finite-time admissibility, without sacrificing both charges, energy, quantum H or equilibrium. Physical pointwise true ηdot(k)>0 must not be turned into a negative right-end value by unconstrained low-degree polynomial projection. Independently validate entire-spectrum and finite-time behaviour, then add coupled χ1/χ2/φ dynamics and realistic SM/FRW. No BEC/relic/CLASS claims.
