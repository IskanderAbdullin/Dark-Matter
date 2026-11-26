# Physical Yukawa 1↔2 invariant event normalization and finite-cell energy gate (2026-11-24 UTC)

**Status:** analytical derivation + independent Python-stdlib scalar moment checks; **not** a constructed nonlinear strong solver, and not an independent numerical 3-species convergence test.

For isotropic, unpolarized, on-shell decay of fermion 1 of mass M into fermion 2 (mass m2) and a real scalar φ (mass mφ), with rest-frame inclusive partial width Γrest (summed final internal states, averaged over initial states), let
`k*=sqrt(lambda(M²,m2²,mφ²))/(2M)`, `Eφ*=(M²+mφ²−m2²)/(2M)`, `E1=sqrt(p²+M²)`, `γ=E1/M`, `β=p/E1`, and `Eφ(p,c)=γ(Eφ*+β k*c)`, `E2=E1−Eφ`, `c∈[-1,1]`. These are **exact** event energies; the phase-space angular distribution is uniform in c only for an isotropic unpolarized rest-frame squared amplitude.

With occupations evaluated at **the actual event energies**, define `F=f1(E1)[1−f2(E2)][1+b(Eφ)]`, `R=[1−f1(E1)]f2(E2)b(Eφ)`. The positive common rate measure per physical volume is

`dW = [g1 Γrest M/(2π² E1)] p² dp (dc/2)`,

and the net event rate is `dRnet=dW(F−R)`. Here g1 counts equally occupied initial internal states and Γrest must be defined consistently for the same summed final channels. No independent g2 or gφ multiplier may be inserted without deriving the channel multiplicity from the matrix element. This normalization follows from `dN1/dt=−g1∫d³p/(2π)³ Γrest M/E1 f1` in vacuum and the isotropic boosted two-body phase space. It is **not** a generic formula for anisotropic/spin-polarized decay or in-medium modified amplitudes.

The weak moments follow from the **same** positive event measure:
`dot N1=−∫dW J`, `dot N2=dot Nφ=+∫dW J`;
`dot ρ1=−∫dW E1 J`, `dot ρ2=+∫dW E2 J`, `dot ρφ=+∫dW Eφ J`.
Thus both number charges and `dot ρ1+dot ρ2+dot ρφ=0` hold pointwise. For `μ1=μ2+μφ` and common T, `F=R` pointwise; in the interior `(F−R)ln(F/R)≥0`. These are **weak continuum** statements, not a proof of an exactly reversible finite-volume strong discretization.

**Independent physical capacities.** For any species s with mass m and occupation per internal state, `V_s(L,U)=g_s[p(U)^3−p(L)^3]/(6π²)` and `W_s(L,U)=g_s[A(U)−A(L)]/(16π²)`, with `p(E)=sqrt(E²−m²)`, `A(E)=E p(E)(2E²−m²)−m⁴ ln((E+p(E))/m)`, and `L>=m`. V has GeV³, W has GeV⁴. For cells crossing threshold clamp L to m. These follow from integrating `g_s E sqrt(E²−m²)/(2π²)` and its E-weighted moment. Numerically near threshold, use stable series/factored differences rather than subtracting nearly equal A values.

**Negative gate:** `W_s/V_s` is a cell-average energy, but for a boosted reaction the actual event energy depends on parent p and c. Depositing number into a single constant-occupation bin and energy as its mean generally fails pointwise `E1=E2+Eφ`; barycentric splitting of daughters can repair weak mean energy while not automatically retaining the same positive inverse quantum factors, equilibrium or discrete entropy. Conversely assigning each event its own independent daughter occupation can violate shared physical DOS and create disconnected modes. Therefore this analytical normalization **does not** establish a positive physical finite-volume strong collision solver.

Independent stdlib check (not CI): `m=86e-9 GeV`, scalar first cell `[m,1e-6]` gives `V=1.6699867896979555e-20 GeV³`, `W=1.2570473805073829e-26 GeV⁴`, `W/V=7.527289366970026e-7 GeV`; cell `[m,1e-4]` gives `V=1.688684520610624e-14`, `W=1.2665138588122812e-18`, `W/V=7.500002773486152e-5`. Reproduce with `python -c` or direct evaluation of the displayed closed forms. No new 3-species rate or CI convergence is claimed.

**NEXT_ACTION:** implement a separately testable positive weak cut-cell integration using this parent normalization, actual Eφ and E2 inside each intersected cell, and exact species V/W capacities. Compare Rφ and all three Q with the existing independently converged analytic thermal reference for pmax=2/8/16 and ≥3 parent/angular/scalar resolutions, including nonzero μ and reverse-dominated cases. Report cutoff and first-cell misses explicitly. Only after this gate attempt reversible physical strong discretization; keep CLASS-ready=false.
