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

## A real 5D Pauli-blocked weak-form quadrature (following event-only audit)

**Pre-registered tests:** (i) correct invariant normalization and both Pauli factors; (ii) common-T zero for every p even with formal nonzero species chemical potentials, (iii) positive entropy and opposite SM heat, (iv) independent rest recoil test, (v) direct strong Cχ(p) integrated number/energy agreeing with symmetric weak moment on **converged finite grids**, (vi) explicit low-cutoff failure saved, (vii) G² and degeneracy scaling. Implementation: [relativistic_weak_form_2to2.py](relativistic_weak_form_2to2.py), [tests](../tests/test_relativistic_weak_form_2to2.py), [reproducible runner](verify_relativistic_weak_form.py), [JSON benchmark](results/relativistic-pauli-weakform-2026-11-21.json).

For isotropic incident χ and e, let P=(E,0,0,p), K=(k,k sqrt(1-μ²),0,kμ), mχ=M, s=M²+2k(E−pμ), v_M=1−pμ/E. A CM outgoing χ is sampled at polar u*=cosθ* and azimuth φ* **relative to the incoming χ CM axis**. It is boosted to lab using (P+K)/sqrt(s). This is not an azimuth-averaged approximation: final e and χ energies depend on both u*,φ*, and their FD/Pauli factors are evaluated at those exact outgoing energies. The spin-averaged tree vector-contact amplitude is the independently established |M|²=2G²[(s−M²)²+(u−M²)²+2M²t], with t=−(s−M²)²(1−u*)/(2s) and dσ/dΩ*=|M|²/(64π²s). The symbol u in the amplitude is the Mandelstam variable, **not** cosθ*.

With F=fχ(E)fe(k)[1−fχ(E')][1−fe(k')], R=fχ(E')fe(k')[1−fχ(E)][1−fe(k)] the *symmetrised* weak energy moment is

```text
Qchi = gchi*ge/2 * integral [d^3p d^3k/(2pi)^6] vM
                      * integral dOmega* (d sigma/dOmega*) (F-R) (E'-E).
Qe = -Qchi; sigma_entropy = Qchi*(1/Tchi-1/Te) >=0
for the FD(Tchi,mu_chi), FD(Te,mu_e) two-temperature trial states.
```

The factor 1/2 is for pairing forward/reverse **weak-form events**. It is NOT a 1/2 identical-particle symmetry factor (χ and e are different). Incoming averaged-spin cross sections are multiplied by gχ ge once. After integrating the trivial initial azimuth and outgoing azimuth normalization, the quadrature coefficient before ∫dp p² dk k² dμ du* ⟨dφ*/2π⟩ is gχ ge/(8π³) * [dσ/dΩ*] * vM * (F−R)*(E'−E). For the *unsymmetrised* strong source at each incoming p, Cχ(p)=ge/(2π)∫dk k²dμ du*⟨dφ*/2π⟩ vM (dσ/dΩ*)(R−F). These are independent formulations: actual finite-grid ∫p² Cχ dp need not conserve χ number exactly without sufficient p/k coverage. We **do not** silently project this residual to zero.

**Measured exploratory benchmark, sandbox Python numerical integration, not BEST or an observational inference.** M=.1 GeV, Tχ=.8 GeV, Te=1 GeV, G=.001 GeV⁻² (nominal contact scale 31.6 GeV if coefficient=1), gχ=ge=4, chemical potentials 0. At n=24 in all five dimensions and pmax=32Tχ, kmax=32Te, Qχ=+2.4742564595×10⁻⁸ GeV⁵ and Qe=−Qχ. A second n=28 grid gives 2.4742562993×10⁻⁸ GeV⁵ (relative change ~6.5×10⁻⁸). This is **numerical refinement only** for the electron-only EFT and isotropic FD trial; it does not calibrate arbitrary M/T=0.1 BGK. It also does not establish full cutoff/UV physical consistency for high-energy tails.

**Necessary negative cutoff result:** a smaller cutoff=15 and n=24 produces weak Qχ=2.450745808×10⁻⁸ but an independently integrated strong Cχ(p) energy moment 2.412972205×10⁻⁸ GeV⁵ (−1.54% relative discrepancy) and a nonzero particle-number source −3.865×10⁻¹¹ GeV⁴. Increasing cutoff to 32 reduces the strong-versus-weak energy discrepancy to ~1.34×10⁻⁷ (relative) and strong number residual to +3.78×10⁻¹⁵ GeV⁴; **a conservative-looking weak integral alone would have hidden this failure**. This is a numerical error due to incomplete phase-space coverage and quadrature; values should be reproducibly regenerated before promotion to verified evidence. At Tχ=Te the paired kernel and Cχ(p) vanish pointwise; direct F−R is separately inspected to catch errors concealed by affinity-based cancellation. Reverse temperature ordering gives Qχ<0 with positive entropy.

Independent normalisation: a **different 2D CM angular/quadrature** for a χ initially at rest (not a thermal gas moment) gives Qrest=1.430815529×10⁻⁸ GeV² at M=1,Te=.1,G=.2; the earlier independent scipy (k,t) integral gives 1.430815515×10⁻⁸ GeV² (~1e−8 relative difference). These are mathematically separate integrations of the same tree amplitude and final electron Pauli.

**Scope gates not passed:** The computation resolves the full quantum collision **weak moment and a sampled strong source only on thermal FD trial distributions**, not a time-evolving source for arbitrary nonthermal fχ(p). Electron mass, other SM species, medium screening, high-T electroweak gauge completion, UV propagator, dark backreaction, realistic gstar(T), simultaneous 1↔2 and Boltzmann-time dynamics are missing. The original phenomenological BGK γ_el is **not replaced**. Thus physical Γ(T), 86-eV ΩDM, BEC, freezeout and CLASS-ready remain unproven. Relevant context: [Binder et al. 2017 kinetic decoupling](https://arxiv.org/abs/1706.07433), [Bringmann/Hofmann 2006](https://arxiv.org/abs/hep-ph/0612238) and [Yoon 2026 BEST](https://arxiv.org/abs/2603.28848), none of which is claimed as independent verification of this particular EFT numerical result.

**Reproduction:**

```bash
python -m pytest -q tests/test_relativistic_weak_form_2to2.py
python kinetics/verify_relativistic_weak_form.py --output outputs/pauli-weakform.json
```

Next scientific gate: replace thermal-FD-only strong source with a generic number-conserving momentum-resolved fχ(q) kernel, construct a positivity- and energy-preserving finite-volume discretisation with explicit SM Q transfer, then compare its linearized relaxation spectrum against BGK. Check also the contact-EFT UV and neutrino/muon channels before cosmological use.

## General nonthermal fχ(p) beyond the two-temperature FD ansatz

The next tested implementation [relativistic_nonlinear_fermion.py](relativistic_nonlinear_fermion.py) evaluates the *same* invariant 5D χe tree collision integral for arbitrary isotropic fermionic χ occupancies, including fχ(p′) in the exact inverse gain and both final χ/e Pauli factors. It now returns actual sampled mode-resolved Cχ(p), its independently integrated number and energy moments, the symmetrised weak heat Qχ=−Qe, and the instantaneous full quantum entropy integrand (F−R) ln(F/R)≥0. No phenomenological BGK projection, no thermal reset. A callable analytic distribution can be supplied. Alternatively `distribution_from_grid(q,f)` uses PCHIP, enforces q=0 and strict monotonic q, preserves 0≤f≤1 on the resolved interval, and explicitly sets f=0 beyond qmax. That tail policy is **a modelling assumption** which is not an actual physical frozen distribution.

An exact analytic FD distribution in this new general kernel reproduces the previous independently implemented FD-only 5D weak Q to roundoff on **identical pmax/kmax**. The common-T FD bath is a pointwise stationary solution, with *all* Cχ(p)≈0 even for nonzero formal chemical potentials. A nonthermal benchmark at M=.1 GeV, Te=1 GeV, G=.001 GeV⁻² and a modulated FD fχ(q)=FD(Eχ/T)[1+ε sin(.7q/GeV)] gives ε=.06, qgrid=321 on [0,34] GeV, and n_p=n_k=24, n_µ=16, n_cm=20, n_φ=16: Qχ≈+5.09292×10⁻⁹ GeV⁵, positive entropy production ≈1.78×10⁻¹⁰ GeV⁴. For small ε the heat responds linearly and entropy quadratically, as required by linear response (reported values are exploratory deterministic quadrature, not a physical relic). Physically charge-conjugate μ prescriptions, thermal masses and SM feedback remain outside this benchmark.

**New negative numerical result:** even when input samples are *exact* thermal FD at Te, PCHIP interpolation and hard high-q truncation can create **spurious heat and entropy**. For 41 tabulation nodes q∈[0,32] GeV, numerical Qχ≈−1.68×10⁻¹² GeV⁵ and entropy≈2.18×10⁻¹² GeV⁴ on a fixed test quadrature. Increasing to 641 points decreases the spurious heat to ~3×10⁻¹⁴ GeV⁵ and entropy to ~3×10⁻¹⁵ GeV⁴ (different grid domains in subsequent scans must be reported explicitly). Exact analytic FD remains strictly stationary. Do not mistake a small nonthermal spectrum source, a tail-induced collision signal or positive numerical entropy for a real thermalization mechanism unless interpolation/conservation error is below the physical signal. With f=0 exactly, numerical log(F/R) is regularized to 10⁻²⁹⁰; the resulting entropy estimator is not a proof of the full H theorem in a truncated phase space.

The operator is **instantaneous** for an imposed fχ(p); there is still no physical FRW time integration or independently calibrated `γ_el` for arbitrary fχ, and no self-consistent SM bath equation. Next scientific work: construct a conservative, quantum-bound-preserving time step using the tabulated mode-resolved collision kernel (without hidden charge correction), check both SM Q and arbitrary-f spectra under q/dt/cutoff convergence; then add gauge-required SM channels and exact in-medium collision physics.

Validation: `python -m pytest -q tests/test_relativistic_nonlinear_fermion.py`. The existing BEST-based FRW pilot is not claimed to have used this new operator.
