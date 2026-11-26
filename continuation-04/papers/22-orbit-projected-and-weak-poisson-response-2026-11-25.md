# 22. Independent orbital susceptibility and regular weak-Poisson response
2026-11-25. Same critical benchmark as 19–21.

## Static result
The fixed-(s,L²) orbit projection P_EL is conditional expectation onto invariants, while P_E projects onto energy alone. Since the latter subspace is contained in the former, the weighted conditional-variance identity gives B_EL[h] ≤ B_E[h] (weight F'(s), nonnegative). Therefore paper 15's analytic bound B_E≤3A/4 also implies A−B_EL≥A/4, including critical radial functions in the same form domain. The static adiabatic response (I−A⁻¹B_EL)⁻¹ has norm ≤4 in the A norm. This does not rule out dynamical resonances or nonlinear evolution, and does not apply to arbitrary nonradial perturbations.

Independent numerical implementation: for s=ψ−q², ℓ=L²=2r²q²(1−μ²), phase volume r²dr q²dq dμ (μ∈[0,1]) equals ds dℓ dr/[8√(ψ−s−ℓ/(2r²))]. For each s find maximum ℓ on circular orbit and both turns for each angular-momentum quadrature point. Transform r=r_p+(r_a−r_p)sin²θ and split at potential FE nodes. Accumulate weighted covariance of potential hat basis on each orbit, then integrate over s,ℓ with F'(s)/Ic. Independently evaluate energy-only covariance and phase measure. Constant outer potential is irrelevant; set h(rt)=0.

Local runs (nf,ns,nℓ,nr per segment): coarse (24,48,16,12), fine (48,96,32,16).
- Linear max B_EL/A .548022 → .552092; energy-shell .553416 → .557703.
- Quadratic .531119 → .544235; energy-shell .536180 → .549919.
- Fine weighted angular-integrated phase-volume discrepancies 8.07e−7 and 2.06e−6; all orbit eigenvalues positive. Min eigenvalue of (B_E−B_EL)/A −2.07e−6 linear and +1.10e−9 quadratic; the small negative is numerical (exact conditional-variance inequality). Further spectral refinement is needed for precision, especially quadratic. Fine inverse gain estimates 2.23/2.19, not certified upper bounds; analytic bound is 4. This independently supports the order-unity radial gravitational feedback, not equality to the finite-time C22 gain.

## Weak causal solver
Replace pointwise enclosed-mass field of paper 21 by FE representation of perturbation U=−Φ, continuous piecewise linear in r², with h(rt)=0 and log-spaced radii. A_ij=∫r²h_i'h_j'dr; exact element coefficient 4(b⁵−a⁵)/[5(b²−a²)²]. Solve A a=Qt∫h δm. The center field is proportional to r in first cell, no ad hoc enclosed-mass projection. A constant shift from exterior shells has no force inside. Independent radial quadrature per element generates the same odd velocity kernel as paper 21; δM(r)=−r² Σ a_i h_i'(r)/Qt maps field coefficients to source. Stream all signed collision tracers and all equilibrium response tracers; causal convolution in field coefficients.

At .5 Gyr, fine n16,nr48,nq32,na16,dt=.001,nf40,r0=.125,ng4:
- Linear Gaussian changes [.0001660265,.0005917469,.0013006733,.0010772724,−.0001143387]; quadratic [.0001823181,.0005617670,.0010291934,.0006450524,−.0001113730].
- Paper 21 pointwise volume fine first four probes differ by ~1–4%; outer small residual more sensitive. Broad amplification persists. Absolute max force .00973/.01203 (rather than spurious O(1–100) values).
- Poisson weak residual 8.94e−17/4.56e−17, total number residual ~1e−18; energy residual −2.08e−7/−1.31e−7. Linear coarse n12,nr32,nq24,na10,dt=.002,nf24,r0=.25,ng4 gives broad [.0001622764,.0005812448,.0012779013,.0010512767,−.0001287179].
- Central force slope linear .03193 coarse → .07918 fine (quadratic fine .10344); relative force outside .1R on fine grid .250/.328. Thus the weak basis removes pathological raw absolute field but does NOT establish converged inner slope or uniform perturbative validity. A finite first cell is a numerical regularization, not physical core modification. Second fine run initially timed out in a combined command and succeeded when run separately.

Reproduce:
python continuation-04/code/lowered_bose_orbit_projection.py --model linear --nf 48 --ns 96 --nl 32 --nr 16
python continuation-04/code/lowered_bose_weak_vp.py --model linear --n 16 --nr 48 --nq 32 --na 16 --dt .001 --nf 40 --r0 .125 --ng 4
Repeat quadratic. Local runs; CI unverified. Next investigate direct pointwise C22 and occupation-relative secular rates in (E,L), especially near-critical core, to prepare evolution rather than extrapolating constant forcing.
