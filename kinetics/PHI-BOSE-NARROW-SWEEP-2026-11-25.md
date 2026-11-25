# Dilute coherent self-resonance: swept band, vacuum particle creation and 4→2 matching
2026-11-25 UTC. Extension of [coherent resonance](PHI-BOSE-COHERENT-RESONANCE-2026-11-24.md) and [scalar completion](PHI-BOSE-SCALAR-QUANTUM-COMPLETION-2026-11-24.md). All results are conditional on an initially homogeneous real scalar, positive independently matched quartic, negligible derivative/higher operators and an unfragmented harmonic action normalized to present Ωh²=.120. Not an actual produced relic.

## New analytical link
For V=m²φ²/2+λφ⁴/24, x=A cn(z|p), p=A²/(1+2A²), the previously derived positive Lamé gap has width Δκ²=3A⁴/4+O(A⁶) around κ²≈3. Independent Duffing expansion in the oscillation phase θ gives x=A cos θ+(A³/16)(cos3θ−cosθ)+O(A⁵), ω²/m²=1+3A²/2+O(A⁴). In the Hill equation, direct cos4θ harmonic coupling of the ±2 Fourier states is +3A⁴/16; eliminating the zero-frequency intermediate yields −9A⁴/16. Net off-diagonal −3A⁴/8 and gap width 3A⁴/4. The center temporal Floquet rate is μc/m=3A⁴/32+O(A⁶). Across the narrow band, μ≈μc sqrt(1-u²), u∈[-1,1]. In adiabatic harmonic expansion A∝a^-3/2, κ∝a^-1, the detuning slope |dκ²/dt|=6H[1+O(A²)]. Thus
I=∫ μdt ≈ (π/4) μc Δκ²/(6H) = (3π/1024)(m/H)A⁸.
The leading approximation requires A≪1 and H≪m; all other interactions, nonlinear depletion and previously produced seeds are omitted.

For linear crossing, canonical two-mode equations dα/dt=−i v t α+g β, dβ/dt=g α+i v t β, I=πg²/(2v). The vacuum mean occupation after one sweep is exp(2I)−1 in the infinite linear-chirp limit; weak production N≈2I. Independent numerical integrations in dimensionless time ±L, with adiabatic positive-norm in/out basis, at I=.01,.1,1 and L=20,40,80 check the formula to <2e−6 at L80; Bogoliubov norm |α|²−|β|² differs from 1 by <1.6e−7. The direct elliptic coefficient monodromy ratios μc/(3mA⁴/32) at A=.03,.05,.1,.2,.3 are 1.00006,.99442,.97791,.91629,.82692; the .03 trace is cancellation-sensitive, not a high-precision proof at smaller A.

## Independent kinetic normalization
Redshifting modes cross k_phys≈√3 m with physical mode flux per volume per time k³H/(2π²). Weak vacuum hard-particle production per volume is therefore (9√3/(1024π)) m⁴ A⁸. Each 4→2 event emits two hard particles, removes four cold particles and reduces total number by two; hard count rate equals net total count loss, not cold loss. With n_c≈6m³A²/λ the incoherent unoccupied-final kinetic event rate is R=√3 λ⁴ n_c⁴/(12288π m⁸), including identical 4!2! factors (Arcadi et al., JHEP 08 (2019) 050, eq 4.15, https://arxiv.org/abs/1906.07659). The weak coherent result equals 2R/4!, exactly. This is the distinct coherent-state local fourth-order correlation, not another symmetry factor to apply to an incoherent Boltzmann gas. The agreement is an independent normalization bridge from a classical background plus quantum vacuum to the perturbative reaction. Neither weak rate can be used as an upper bound on stimulated production at I≫1. A separately matched local sextic changes the amplitude; the pure quartic approximation here is consistent with neglect of the tiny previous finite fermionic c6 contribution for these two λ examples.

## Conditional crossing epochs
Radiation H from photons, massive e± Fermi integrals and 3 instantaneous-decoupling neutrinos (Tdec=2 MeV, no QED/reheating correction), SM entropy sets a(T); scalar normalization Ωh²=.120. The diagnostic assumes uninterrupted homogeneous background, although earlier resonance may already fragment it. T is visible photon temperature; q=√3ma is today's comoving momentum (a0=1).
| λ | I | T (keV) | A | q (eV) | incoherent net Γ/H | vacuum hard-on-cold elastic Γ/H |
|---|---:|---:|---:|---:|---:|---:|
|4.698884e−13|.01|215.809|.004657|1.28458e−7|4.56e−10|.431|
|4.698884e−13|1|323.938|.009268|8.11944e−8|1.15e−8|.694|
|4.698884e−13|10|400.746|.013078|6.45353e−8|5.79e−8|.879|
|3.82587e−8|.01|2.97294|.001516|1.17677e−5|3.51e−4|2.95e7|
|3.82587e−8|1|4.71180|.003024|7.42490e−6|.00881|4.68e7|
|3.82587e−8|10|5.93180|.004272|5.89781e−6|.04416|5.89e7|

The vacuum hard-on-cold elastic diagnostic uses σ=λ²/(32πs), s=6m² for E_h=2m, v=√3/2, Γ=n_cσv. Its ratio to incoherent net number loss is 16/(λ n_c/m³)². Hard elastic scattering can be very fast on an expansion timescale even when the perturbative number-changing rate is not. It does not establish a kinetic relaxation coefficient in a coherent medium. The large-λ branch retains formal swept resonance even when the emitted q exceeds the previously selected 1.14255e−6 eV rms reference; this does not prove a final hot fraction. The small-λ example's I=1 q is substantially below that reference. No estimate of final PSD, Ω, BBN likelihood, or exclusion. The underlying coherent background is conditional, not a production mechanism.

Reproduce: python kinetics/phi_bose_narrow_sweep.py. Local SciPy checks: six crossing roots, five exact Lamé monodromies, nine linear chirp solutions, independent analytic normalization; all stated checks passed in local independent runs. New remote CI not checked. Next: distinguish residence-time scattering from post-sweep redistribution and estimate when formal linear emission invalidates the background; then nonlinear kinetic/field treatment or independent halo constraints. CLASS-ready=false.
