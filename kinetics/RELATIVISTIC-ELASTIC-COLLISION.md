# Relativistic χ–electron finite-temperature collision: invariant formulation and event audit

**Status 2026-11-21:** first **on-shell event-only** numerical gate, not a solved 2→2 collision integral, not a physical gamma_el and not a validated 86-eV relic. Natural units ħ=c=kB=1. This document explicitly prevents promoting the old heavy-NR drag to M/T=0.1.

## Hypothesis and falsification gates (fixed before event audit)

For Dirac χ and e with vector contact G(χbar γμ χ)(ebar γμ e), a correctly normalized Lorentz-invariant collision integral should conserve elastic χ/e numbers and total four-momentum, vanish at common T with arbitrary species chemical potentials, and produce nonnegative entropy. Negative control: Tχ≠Te generally gives nonzero local gain/loss. The contact approximation requires momenta and temperature well below the mediator scale; the electroweak symmetric phase requires gauge-completed e_L,ν_L,e_R operators, with separately derived couplings and degeneracies.

## Complete formal collision source (not yet integrated)

Let dΠ_i=d³p_i/[(2π)³ 2E_i], P+K=P′+K′ and |M|² be spin-averaged with spin/antiparticle degeneracies accounted consistently. For one χ spin state (matrix-element spin conventions must be adjusted accordingly),

Cχ(p)= [1/(2E_p)] Σ_{e−,e+} ∫dΠ_k dΠ_p′ dΠ_k′ (2π)⁴δ⁴(P+K−P′−K′) |M|²
× { fχ(p′)fe(k′)[1−fχ(p)][1−fe(k)] − fχ(p)fe(k)[1−fχ(p′)][1−fe(k′)] }.

The corresponding electron operator interchanges χ and e labels, using the **same** transition kernel. This guarantees ∫d³p/(2π)³ Cχ=0, ∫d³k/(2π)³ Ce=0, and Qχ+Qe=0 with Qχ=∫d³p/(2π)³ Eχ Cχ, provided C denotes the rate contribution to ∂t f (otherwise account for the usual 1/E Boltzmann streaming convention). The collision measure and internal-state averaging must not be double-counted. For each event F=fχ fe(1−fχ′)(1−fe′), R=fχ′fe′(1−fχ)(1−fe), and (F−R) ln(F/R)≥0. With FD populations at common T, F=R for any fixed μχ, μe because elastic reactions preserve both numbers.

**Caution about normalizations:** the invariant expression above with prefactor 1/(2E) is the collision contribution to ∂t f for a homogeneous physical-momentum distribution. The prior published |M|²_avg averages incoming spins; to compute a single-spin f source one must sum over electron initial spins and final spins while holding the χ initial spin fixed, or consistently use unpolarized species distributions and matching degeneracy factors. These are not yet implemented as a full kernel.

For massless incoming e, s=M²+2k(Eχ−p cosθ), CM momentum p*= (s−M²)/(2√s). Construct outgoing P′* with energy √(M²+p*²) and isotropic angles relative to incoming χ*; boost with β_CM=(pχ+k)/ (Eχ+k), and K′=P+K−P′. The physical differential weight is dσ/dΩ*=|M|²/(64π²s) for elastic initial/final CM momentum equality, with unpolarized |M|² from MICROPHYSICS-MODEL.md; integrate k²dk dcosθ dΩ* and all final FD factors, **without** replacing them by an azimuth-averaged unblocked energy kick. The moving χ gain term additionally needs inverse transitions and phase-space Jacobian, not merely a test-particle loss rate.

## Reproducible numerical event gate

Command: `python kinetics/relativistic_elastic_audit.py`. Deterministic Python stdlib random seed 21112026, 4000 nonuniformly sampled kinematic events, M and T each log-uniform in [0.1,10] arbitrary common energy units, p/T,k/T in [.05,5], incoming angle in [.02,3.12], outgoing CM polar angle in [.02,3.12], random μχ,μe/T in [-.2,.2], Tχ/Te in [.7,1.3]. Four-vectors are boosted directly, and gain/loss computed separately from FD factors. The implementation does not use the earlier event ledger code and tests a different kinematics construction.

Locally executed in this research run: maximum relative on-shell error 3.6938160863364544e-13; four-momentum residual 0; common-T relative gain/loss error 1.3319660326089454e-14; minimum sampled nonequilibrium event entropy 4.323534507669499e-50; weighted heat/entropy-sign sum 0.42937267568510806 (not an integral; no physical units); 4000 samples. Note these values refer to the independently run pre-commit local version of the same algorithm; rerun the committed script and record exact output/CI before using as an independent quantitative benchmark. Checks use numerical floating-point thresholds 1e-12,1e-11,1e-15.

**Negative result:** This test does not establish a physical thermalization rate, equilibrium of the *integrated* operator, finite-temperature plasma dynamics, or the relativistic BEST benchmark. It checks the necessary event geometry and detailed-balance algebra only. Random event sampling is not a converged quadrature, and no collision-kernel discretization has yet passed an H theorem.

## Immediate next experiment

Implement a gain/loss weak-form collision integral in which each forward event contributes Δψ=ψ(P′)−ψ(P) to the χ moment with full Pauli F−R, using the exact CM differential weight and initial χ,e phase-space measures. Derive and test all 2π and degeneracy factors independently against the rest unblocked rate; then run deterministic Gauss-Legendre/angular refinement (including outgoing azimuth) and two independent quadratures. Test ψ=1, E, pz and electron counterparts, common-T null, temperature-reversal Q, and entropy. Only after this may one project a linearized kinetic relaxation eigenmode; it need not be describable by a single BGK gamma.
