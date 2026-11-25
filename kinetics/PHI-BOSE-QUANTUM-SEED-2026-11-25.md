# Resonant-band vacuum normalization and nonlinear fragmentation (2026-11-25)

## Physical normalization, explicitly conditional initial state

At a turning point of the prescribed homogeneous field, define instantaneous positive-frequency oscillators with dimensionless frequency ω_k=√(k²+a+3b). With z=m√(1+2A²)t, physical background amplitude φ_max=√12 m A/√λ, C=[m√(1+2A²)/φ_max]²=λ/(6b), equal-time symmetrized covariances per mode are C/(2ω_k), C ω_k/2, vanishing cross covariance, with ∫d³k/(2π)³ measure. Evolve each two-dimensional covariance via exact linear monodromy M V Mᵀ. Count only the finite resonance band; subtract its initial quadratic zero-point energy, normalize to E_bg=a/2+b/4. No UV vacuum populated; this is NOT a renormalized interacting vacuum, full quantum state or relic prediction. 40/80 Gauss-Legendre nodes agree at the reported precision; max |det M−1|<5.2e−10.

| λ | A at illustrative T=20 MeV | band seed rms | initial band energy/E_bg | first linear ΔE/E_bg>1% (periods) | first >10% | H Δt at >10% |
|---:|---:|---:|---:|---:|---:|---:|
|4.698884e−13|2.72085|1.14017e−8|2.23487e−15|70|76|2.90e−13|
|3.82587e−8|120|3.32938e−6|2.04703e−10|37|42|3.81e−15|

Radiation SM toy background from prior sweep; no thermal backreaction. At the linear 10% crossings, z=555.762 and 311.482, sigma=.08491 and .08308; linear approximation already marginal. Expansion is negligible *during this local instability* if these coherent initial conditions occur.

## Classical stochastic follow-up

In periodic side L=32 box, sample independent real/imaginary Gaussian coefficients in each independent Hermitian pair with variances C/(4ω L³) and Cω/(4L³). Seed=42. Same realization on both grids. Use dealiased nonlinear solver with 128 steps/period; finite resonant band only, so approximation is classical-statistical in highly occupied amplified modes, not full quantum scattering. Initial sampled sigma from direct band variance 1.08975e−8 and 3.00374e−6; fixed cancellation bug in sigma diagnostic for fluctuations below floating point resolution by directly summing nonzero Fourier modes (def6b10). This does not affect evolution.

| A | N | periods | first D>10% | max sampled D | final D | final U(k>2) | max sampled energy error |
|---:|---:|---:|---:|---:|---:|---:|---:|
|2.72085|48|105|77|.31606|.22825|.24846|.0008683|
|120|48|70|42|.54884|.45914|.36465|.0009400|
|120|64|70|42|.53345|.45243|.40170|.0009514|

D=1−E_mean/E; U is high-k share of nonzero-mode kinetic+gradient energy, not particle fractions. Large-lambda nonlinear onset agrees with continuum linear estimate; smaller-lambda differs by one period. D subsequently oscillates; do NOT interpret as irreversible depletion or stable condensate fraction. Spectrum not cutoff-converged. At T=20 MeV conversion q/k_dimensionless=2.86588e−9 and 1.22335e−7 eV; this is kinematic conversion only, not final q_rms. No fermions, expanding background, full quantum collision, physical relic, or observational exclusion. In particular, fragmentation invalidates uninterrupted homogeneous normalized action trajectory when such coherent initial conditions are imposed; does not demonstrate an initial field existed.

Reproduce from kinetics: python phi_bose_quantum_seed.py --linear ; python phi_bose_quantum_seed.py --N 48 --periods 105 --save quantum_mixed48.npz ; repeat --A 120 --lam 3.82587e-8 --N 48/64 --periods 70 (distinct checkpoint paths). Source e3b6121. Local checks only; remote CI not checked.
