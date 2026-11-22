# Backlog — prioritize physical bottlenecks, not commit count

This is a living scientific queue. For the next executable task see [NEXT_ACTION.json](NEXT_ACTION.json). After every completed step, recompute priority by (physical impact × reduction of uncertainty × independent testability) / cost; do not equate `CI=green` with `physics=validated`.

## P0 — physics necessary before credible 86-eV freezeout

- [ ] Choose physically coherent χ1/χ2/φ Lagrangian and SM portal; derive `Γ(T)`, thermal masses, elastic `γ_el(p,T)`, Q, detailed balance, perturbativity and stability.
- [ ] Quantify importance of other stimulated 2→2 and SM number-changing channels, justify inclusion/neglect for each parameter regime.
- [ ] Replace phenomenological BGK relaxation with physically calibrated operator and independent moment/source benchmarks.
- [ ] Couple `g*ρ(T)`, `g*S(T)`, SM heat/entropy, Friedmann H and all three species including IR zero mode.

## P1 — trustworthy late-time spectra

- [ ] Integrate from thermal initial conditions to genuine relic freezeout/saturation, not a fixed early `tΓ`.
- [ ] Joint convergence in `dt`, q-grid, `qmin`, `qmax`, MC seeds/neval; error-budget and independent modes/energy/Q checks.
- [ ] Determine whether `n0` remains nonzero as `qmin→0`, separately from `N0∝qmin³`; quantify hot/warm component.
- [ ] Compute absolute dimensionful Ωφ and distribution without renormalization to ΩDM.

## P2 — observational predictions and independent publication-grade audit

- [ ] Validated transfer/perturbation calculation (CLASS or justified alternative), physical treatment of zero-mode/residual collisions.
- [ ] Compare CMB, BBN/ΔNeff, Lyman-α, LSS, relevant decay/interaction limits with source data and systematics.
- [ ] Independent rerun/model comparison; final allowed/excluded parameter region or rigorous inconclusive constraints; reproducibility package.

## Already supported but **not** final

- [x] Actual pinned external BEST 1↔2 short benchmark and finite massive φ diagnostics.
- [x] Corrected short BEST+SM BGK finite-heat conservation with identical quadrature weights ([run](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44929647245)).
- [x] Preserve negative old quadrature plateau; do not treat early `N0∝qmin³` as BEC.

## Execution conventions

Every stage: hypothesis → control/independent proof → real measured experiment → uncertainty and negative test → checkpoint → commit/push → immediately choose the next stage. If one item is blocked, select another productive item and record blocker in STATE/NEXT_ACTION. Do not use PC/Computer Use without explicit user approval. Heavy public Actions only when justified; keep artifacts small.

## New completed A-phase substeps (21 Nov 2026; not a completed P0)

- [x] Derive explicit tree-level scalar-Yukawa vacuum width, threshold dependence and BEST sum-normalization crosscheck ([model](../MICROPHYSICS-MODEL.md), [CI #44931903023](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44931903023)).
- [x] Derive vector EFT massless-electron total and transport cross sections, independently integrate dσ/dt, rest electron transport collision frequency and heavy-target asymptote; explicitly reject it as gamma_el.
- [x] Derive on-shell FD/BE quantum detailed-balance identity for decay vs inverse, dimensional BEST coupling bridge ([CI #44932046451](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44932046451)).
- [x] Implement χ-at-rest Pauli-blocked electron recoil heating and independent angular/heavy χ controls (CI #44932154301 was running when checkpoint written).
- [ ] **Immediate next:** complete moving χ thermal elastic operator, gauge-required neutrinos, final-state χ blocking and full SM transfer, then independently calibrate kinetic relaxation/drag gamma vs BGK. Do not substitute pilot observables.

These completed substeps are tree-level partial models. Full first P0 item (UV-consistent portal, thermal masses and physical rates) remains unchecked. Check CI before promoting E-009.

## New completed A-phase steps (continued, 21 Nov 2026)

- [x] On-shell elastic χe event gain/loss with both Fermi blocking factors, common-temperature detailed balance, entropy and energy exchange sign; no full physical integral yet ([CI](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44932341208)).
- [x] Moving-χ incoming flux and transport collision frequency ν_T(p,T), rest p→0 verified; still **not** a thermal gamma ([CI](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44932437443)).
- [x] Exact outgoing CM angular moment: moving χ signed laboratory energy heat Q(p) without final Pauli; positive for slow and negative for fast χ, independent t-integral validation ([CI](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44932561256)).
- [x] Nonrelativistic *conditional* Fokker–Planck eta_NR=Q_rest/(3T) from microscopic Pauli-blocked recoil and Einstein fluctuation–dissipation, with hard rejection for M/T<100; CI #44932655044 completed successfully. **Not applicable to historic M/T=.1**.
- [ ] **Priority now:** full relativistic, final-state-Pauli electron+neutrino collision operator, physical matrix element, gauge-qualified UV validity and moving χ momentum/heat relaxation; only then couple to BEST/SM without phenomenological gamma.

- [x] Relativistic on-shell event geometry and FD detailed-balance audit (4000 seeded local cases), plus formal invariant 2→2 gain/loss; **NOT** an integrated operator ([report](../RELATIVISTIC-ELASTIC-COLLISION.md)).
- [ ] Next: integrate invariant weak-form Pauli collision moments with actual vector-contact |M|², independent angular refinement and all normalization/degeneracy factors.

## 2026-11-21 next 5D quantum thermal collision milestone

- [x] **Full final-state Pauli 5D invariant weak 2→2 χe collision moments for isotropic FD(Tχ), FD(Te), and sampled unsymmetrised Cχ(p)**: temperature reversal, common-T fixed point, entropy sign, explicit Møller flux, CM outgoing azimuth and spin-degeneracy convention, and independent rest recoils. [CI success #44935320673](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44935320673); [code](../relativistic_weak_form_2to2.py); [JSON](../results/relativistic-pauli-weakform-2026-11-21.json).
- [x] **Negative quantitative discretisation control**: at cutoff15 strong Cχ number drifts and independent strong/weak energy differs by −1.54%; enlarging to cutoff32 reduces discrepancy to ~1e−7. Do not declare conservation from weak-form algebra alone.
- [ ] **Immediate next**: extend from FD-only trial to arbitrary nonthermal χ occupations; explicit q-grid interpolation, Pauli fχ′, sampled Cχ(p), independent number/heat/entropy tests, and a genuinely conservative finite-volume time step. Then decide if and when a BGK approximation is defensible.
- [ ] Electron-vector EFT still incomplete: full gauge/UV completion, muon/neutrino channels at physical T, finite me, medium effects, thermal 1↔2 collision, realistic gstar/H, and 86-eV late relic remain unvalidated.

## 2026-11-21 full-Pauli nonlinear distribution follow-up

- [x] **Replace FD-only assumption with actual arbitrary isotropic fχ(p)** in invariant 5D e± collision, including incoming/outgoing χ and electron Pauli, reverse processes, strong Cχ(p), symmetric weak Q and instantaneous entropy rate. Analytic FD comparison passes locally (code [here](../relativistic_nonlinear_fermion.py), CI #44935768129 pending in checkpoint).
- [x] **Quantify false heat due to thermal FD interpolation**: PCHIP81 at qmax34 produces spurious Q~−2.57e−12 GeV5 and entropy~6.55e−14; PCHIP641 residual Q~+6.62e−15, analytic FD ~1.39e−24. Preserve, do not turn these into physical thermalization signatures.
- [x] Perturbed FD with nonthermal ε=.03/.06/.12 yields approximately linear heat and quadratic total-entropy production as ε→0. 5D quadrature and qgrid refinements recorded in [JSON](../results/relativistic-nonlinear-pauli-2026-11-21.json).
- [ ] **Next**: genuine time-evolving, Pauli-bounded, number- and energy-conservative finite-volume collision network; 2-sector SM Q and independent weak-form tests. SM UV/gauge/thermal rates and cosmological goals still unmet.

## 2026-11-21 next corrections and thermal moment benchmark

- [x] FD logit-energy interpolation eliminates numerically spurious PCHIP equilibrium Q even on very coarse q grid; conservatively deposit **moments** of microphysical events to neighboring energy centers with exact discrete charge/energy ledger. **Not a Pauli-stable time integrator** ([CI #44936130211](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44936130211)).
- [x] Minimal SM SU(2)_L lepton completion pure electron vector χ contact implies neutrino channel. For a *common massless bath*, G_L=G_R=G gives χ collision weight ×6/4=1.5 vs electrons only. This is not a UV-complete plasma model ([CI #44936289831](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44936289831)).
- [x] Number-conserving FD(Tχ,μ(Tχ)) thermodynamic heat capacity and conditional **thermal energy-moment response** from physically integrated 5D full-Pauli kernel; gamma_T,moment≈3.03183e−7GeV for M=.1,T=1,G=.001. CI #44936484677 pending at checkpoint. It is not the true kinetic eigenmode/physical γ_el.
- [ ] **Now:** compute and independently converge linearized nonthermal Pauli collision operator/eigenspectrum, compare heat mode vs fixed-N γ_T. Then build actual Pauli/charge/energy-preserving finite-volume time stepping and coupled finite SM Q ledger. Continue to realistic UV/thermal channels and eV freezeout only if validated.

## 2026-11-21 linearized quantum relaxation basis audit

- [x] Reproduce one fixed-N FD thermal energy moment from **independent** 5D Pauli Dirichlet two-mode matrix, susceptibility, common-T stationary point and centered ±ε nonlinear collision source (CI [#44936880435](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44936880435)). A restricted thermal moment is NOT universal relaxation gamma.
- [x] Independent orthonormal FD logit polynomial basis 2→5 modes and n20/24 full 5D collision quadrature demonstrate **negative basis-convergence result**: smallest Ritz 2.125e−7→1.329e−7 GeV (−37.45%), whereas collision n20→24 ≈0.019%. [Model](../higher_order_pauli_modes.py), [compact scan](../results/pauli-ritz-basis-scan-2026-11-21.json); CI for commit `bd776da0` must be checked.
- [ ] **Next:** refine true slow mode using 6+ polynomial and localized IR/UV perturbation basis and/or variational bounds. Separately do real Pauli-bounded finite-volume time dynamics plus SM heat. No physical 86eV relic/CLASS until realistic UV, thermal rates and long-term freezeout.

## 2026-11-21 finite-SM quantum time-integration safety benchmark

- [x] Derived bounded scalar implicit solution for each pairwise Fermi transition with ΔN bracket ensuring 0≤f≤1, finite ESM=A T⁴, QSM=−Qdark, χ number, exact energy and total entropy monotonicity; 130 random cases, endpoints, exact FD fixed point, stiff dt 1e5. [Code](../fermion_finite_bath_pairs.py), [report](../FINITE-BATH-QUANTUM-PAIR-REFERENCE.md), CI for `c373bca1` pending.
- [x] Negative time-accuracy control: split 3-bin network passes all invariant/entropy tests but 10 vs 100 substeps at dt=1 differs in f by >0.001; time convergence remains necessary.
- [ ] **Physics blocker P0 still open**: link reduced numerical transition Γ_ij to the **actual** full-Pauli 5D χe scattering kernel and gauge-required νe bath, with electron recoil/chemical µ and Pauli, physical phase-space cells; avoid identifying toy Γ with a derived kinetic relaxation gamma.
- [ ] Then true finite-SM Q, gstar(T), Friedmann H, actual BEST Bose 1↔2 and late 86eV relic/zero-mode/observations. CLASS-ready=false.

## 2026-11-21 further collision spectral audit and next physics

- [x] Independent **three-term Stieltjes** recurrence avoids unstable high-degree monomials and matches earlier 2/5-mode Pauli matrix eigenvalues at same 5D mesh. [Test](../../tests/test_stieltjes_pauli_modes.py); CI must be checked.
- [x] New scientifically important **negative result**: smallest nested full-Pauli Ritz rate falls through 11 modes: 2.12537e−7 (2) → 1.32934e−7 (5) → 8.56476e−8 GeV (11), while n20/24 grid changes only 0.041% for 11 modes. No demonstrated physical kinetic spectral gap or reliable scalar gamma yet ([report](../PAULI-STIELTJES-11MODE-AUDIT.md)).
- [ ] **Next priority**: local IR/high-energy packets and rigorous variational/tail bounds; convert 5D χe Pauli transitions into physically derived finite-volume conserved Fermi capacities, then finite SM heat and real time evolution. Do not count arbitrary Γ in reduced network as microphysical progress.
- [ ] True cosmological 86-eV freezeout/Ω/BEC/CLASS and peer-independent observational test remain full future stages.

## 2026-11-21 localized IR variational collision experiment

- [x] Explicit susceptibility-number-orthogonal Gaussian χ logit packet and full 5D Pauli collision Dirichlet Rayleigh with exact CM/electron recoil. [Code](../localized_pauli_packets.py), [report](../PAULI-LOCALIZED-PACKET-AUDIT.md).
- [x] **Negative quadrature discovery**: uniform 20-node χ p-grid spanning [0,32T] falsely underestimates narrow IR Rayleigh by ~46% because it fails to resolve p/T≈.23; split p/k at 2T stabilizes n20→24 to ~0.069% and cutoff20→40 under 0.004%.
- [x] IR packet E/T=.25,width=.08 R≈4.85e−8GeV below 11-mode polynomial Ritz 8.56e−8GeV, indicating still insufficient global polynomial basis. This is NOT the true kinetic eigenmode and no full spectral gap proof.
- [ ] **Physics P0**: derive actual microphysical 5D Pauli χe transitions into finite-volume Fermi capacities and finite SM heat, compare to independent weak/strong, implement time integration. Constant toy pair Γ remains disallowed as physical γ_el. Gauge/UV cosmology and 86-eV Bose relic still open.

## 2026-11-21 finite-bin local detailed balance gate

- [x] Derive per-state Kij/Kji thermal ratio, FD fixed point, finite-bath energy and entropy pair identities; independent seeded synthetic 10k algebra control (not χe physics). [Report](../BIN-LOCAL-DETAILED-BALANCE.md).
- [ ] Integrate **actual** paired 5D χe rates on adapted IR physical phase-space bins, verify local detailed balance vs bin-center projection and independently converge heat/strong Cχ; if impossible with fixed center energies, retain negative result and use energy-resolved conservative discretization.
- [ ] Pauli-bounded physical finite-SM evolution and realistic gauge/UV/thermal channels; no 86-eV relic or CLASS claim.

## 2026-11-21 paired finite-bin recoil counterexample

- [x] Synthetic continuous event measure with exact pointwise detailed balance exhibits converged finite-cell center-energy balance bias O(h²), +2.277% at width .8; [report](../FINITE-BIN-RECOIL-PROJECTION-CONTROL.md), [script](../finite_bin_recoil_projection_control.py). **Not physical χe 5D rates.**
- [ ] P0 unchanged: evaluate actual paired 5D χe recoil event integrals with physical Pauli and phase-space capacities, quantify correction and heat error, then finite-SM time integration.

## 2026-11-21 full 5D χe microphysical finite-bin experiment

- [x] First real physical 5D χe Pauli tree-level scattering event graph with exact radial p³-state capacities and positive outgoing energy-barycentric deposition. Physically normalized equilibrium energy-mode D agrees with independent continuum 5D source within 0.0386% at χ128 (specific EFT benchmark); see MICROPHYSICAL-BIN-GRAPH-AUDIT.md.
- [x] **Negative result:** event-energy mean conserving does not preserve quadratic Dirichlet: +7.18% at χ16, +0.41% at χ64, +0.10% at χ128 from strictly positive Jensen variance; asymmetry of physical raw L creates spurious common-T FD heat even when true continuum Q=0.
- [x] Symmetric finite-bin L graph can be used as a **declared linear-response model closure** with positive local Pauli pair rates at reference T, but it is not a fully paired nonlinear electron/χ Boltzmann kernel or physical all-T Γ.
- [ ] **P0:** construct paired energy-conserving χ+e collision quartets with actual microphysical transition weight and full Pauli, no posthoc symmetry; compare physical Q and entropy/capacity and evolve with finite SM bath. Then EW-gauge/UV/Bose/FRW relic/CLASS only on validated source.

## 2026-11-21 real χ+e energy-conserving 5D quantum collision quartet

- [x] **Actual physical amplitude 5D χe events → shared χ/e energy lattice** with equal center total energy, physical radial capacities, complete χ and e Pauli forward/reverse, separate number and energy ledgers. Finite-bin rate tensor Wraw is not symmetric; Wpair=(Wraw+Wreverse)/2 is a declared model-changing approximation, never mislabel as exact 5D continuum. [Code](../paired_energy_quartets.py), [report](../PAULI-CHI-E-QUARTET-AUDIT.md).
- [x] Independent nonthermal Gauss 5D weak heat reference: for χ36/e54/dE=.5/nangular10, Qχ(.9)=+2.0155e−8 vs +1.9853e−8 and Qχ(1.1)=−4.4821e−8 vs −4.4325e−8 GeV5 (~1–2% deviations). FD-weighted omitted collision rate ~4.16e−5; grid accuracy NOT universally converged.
- [x] **Negative:** raw microscopic finite energy-lattice W Pauli-weighted forward/reverse asymmetry 3.69%, false thermal Q+5.795e−12 GeV5; averaging W fixes exact FD but changes physical rates. Preserve and use as primary falsification gate.
- [x] A single true 5D-derived *pair-averaged* quartet supports implicit physically bounded χ+e fermion time transfer conserving species numbers and total energy with entropy growth; this is NOT a full multi-quartet evolving χe bath.
- [ ] **Immediate P0:** directly time-reversal-paired 5D finite-cell quadrature without posthoc W averaging; sparse all-quartet Pauli-safe quantum time integration with measured dt/grid convergence; full SM ν/e plasma and BEST Bose 1↔2 remain open.

## 2026-11-21 full finite χ/e Pauli quartet network time evolution

- [x] Actual coupled χ/e occupation evolution through all 1092 positive-rate energy-changing quartets for small physical 5D-derived (but explicitly rate-pair-averaged) lattice. No fixed electron thermostat. Each implicit Pauli substep preserves χ/e numbers, total energy and entropy monotonicity; globally verified at dt2e6GeV^-1 for 1,4,16,64 substeps. [Source](../paired_quartet_network.py), [time audit](../PAULI-QUARTET-NETWORK-TIME-AUDIT.md).
- [x] Independent initial instantaneous Pauli Qχ and small-dt numerical derivative agree: dt1000 relative error 0.0347%; large dt nonlinear changes expected. Thermal common-T independent μχ/μe fixed point holds.
- [x] **Negative:** exact conservation and positive entropy do NOT establish time/order convergence. One-step ordering discrepancy 1.44e−4 GeV4, shrinking to 1.35e−6 after 64 substeps; underlying raw W still requires a rate-changing symmetrization. No physical Γ_el, Bose freezeout or CLASS.
- [ ] Immediate physical priority: derive directly microreversible paired 5D on-shell event quadrature without ad hoc W averaging, independent nonlinear χe q/t/cutoff/angular convergence and realistic SM lepton heat. Then pinned BEST Bose 1↔2 and 86-eV cosmology.

## 2026-11-21 physically direct microreversible χe q-transfer kernel — no W averaging

- [x] Derive on-shell four-body Lorentz invariant χe 2→2 matrix-element phase-space angular convolution as **positive one-dimensional momentum-transfer integral**, fully analytic cone azimuth average and correct (gχge/(256π^5)) spin-normalised prefactor; directly reciprocal under event time reversal, no Wraw→(Wraw+Wraw^T)/2. [Report](../INVARIANT-CHI-E-ENERGY-KERNEL.md), [source](../microreversible_chie_lattice.py).
- [x] Populate χ+e shared energy lattice from physical exact K(E,k;E',k'), retaining exact center energy and correct state capacities; nonlinear full-Pauli source and implicit all-quartet time stepping conserve χ/e charge, energy and entropy with proper Pauli bounds. **Energy collocation still a discretization**, not full continuum kinetic evolution.
- [x] Independent 5D Gauss CM n28 heat vs direct ℓ integral plus energy lattice: h=.5, cold/hot relative discrepancies ~0.017% each, h=.125 <0.007%; massless-e conditional benchmark only. Negative h nonmonotonicity preserved: h=.25 errors ~0.08–0.09% exceed h=.5 errors.
- [ ] **Verify GitHub CI** after Python syntax fixes: old f172 and 916f failures were engineering regressions, corrections are a3aa and bcd47b; do not promote until all new tests pass and runner reproduces compact data.
- [ ] Physics P0: finite electron mass and e± thresholds/annihilation, electroweak lepton doublet ν channel and UV mediator validity, real gstar and finite SM heat; couple genuine microphysical source to pinned BEST 1↔2 Bose/Pauli and FRW only after independently validated plasma/bath and h/dt/cutoff convergence. No 86-eV freezeout, ΩDM or CLASS claimed.

## 2026-11-21 finite electron mass and CI gate

- [x] Verified actual earlier [Scientific CI #44947528381](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44947528381) success for HEAD 078c3c4 (previous heat float tolerance correction).
- [x] Derived and committed conditional finite-mass electron invariant ℓ kernel and massless/reversal/Pauli/forbidden-transfer regression; local standalone controls passed, **new repository CI still pending**. [Derivation](../FINITE-MASS-CHI-E-KERNEL.md).
- [ ] Verify new real CI, independently derive massive spin trace, compare full massive 5D CM heat vs ℓ kernel at T/me=0.1,1,10 with at least three resolutions and electron cutoff; preserve negative results.
- [ ] Build finite e± plasma with number-changing/chemical potentials and finite heat capacity, gauge-consistent portal and EFT temperature validity before cosmological evolution.

## 2026-11-22 independent finite-electron-mass 5D CM and spin-trace audit

- [x] Previous finite-me [Scientific verification #44954815892](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44954815892) **success** on 065644e; new benchmark commit CI not yet certified.
- [x] Compute **independent** CM 5D p,k,μ,cosθ*,φ* with actual Lorentz boosts and Fermi forward/reverse, plus separate exact ℓ-cone energy integral at Te/me=.1,1,10, Tχ/Te=.9/1.1. Fixed n40/n160 comparison <8.12e−7 relative in all six benchmark cases; test hot and cold sign. [Audit](../FINITE-MASS-CHI-E-5D-INDEPENDENT-AUDIT.md).
- [x] Explicit independent Dirac 4×4 trace (four unequal-mass physical kinematics) matches finite-mass analytic |M|² to ≤1.5e−15 relative; massless previous 5D normalization reproduces.
- [x] **Negative:** E/ω/Δ quadrature convergence nonmonotone and cutoff sensitivity at ~1e−5; n160 coincidence is not a uniform error bound. Full physical SM source and thermal rates remain unknown.
- [ ] **Next P0:** verify CI for new 5D tests; tighten cutoff/crossing with e−/e+ chemistry and charge neutrality, finite γ/ν SM heat and e+e−↔γγ, plus explicit gauge/UV completion and actual G(T,q) validity. Follow with nonthermal Pauli χ spectra and time convergence coupled to physical finite bath. No physical gamma_el/BEC/86eV relic/ΩDM/CLASS.

## 2026-11-22 finite e±/γ thermal response and neutrality

- [x] Derive finite-mass ideal e−+e++γ common-T equilibrium EOS with ρ, P, s, C at μ=0 from explicit FD phase-space and check against independent fugacity Bessel K₁/K₂ expansion and fixed-m finite-difference derivative. [Audit](../EM-PAIR-THERMODYNAMICS-AUDIT.md), [code](../em_pair_thermodynamics.py).
- [x] **Negative:** massless g*ρ,EM=5.5 at T/me=.1 is wrong (actual 2.015461); charge-neutral illustrative η=6e−10 at T/me=.025 needs μe−/T=13.8377 and suppresses e+ by 10⁻¹² relative to e−. Explicitly NOT cosmological η evolution or an actual annihilation collision calculation.
- [ ] Verify CI for both new massive 5D and e± EOS reports and independent regressions; 5D run #44956695501 is now fully successful on Python 3.11 and 3.13; next EM EOS CI pending.
- [ ] **P0:** finite charge-constrained C and e+e−↔γγ Pauli/Bose matrix-element/rate, photon entropy/energy and neutrino weak-sector ledger, evolving SM heat without bath reset; gauge-consistent χ electron/ν portal and full nonthermal χe kinetic eigenvalue, then actual Bose decay FRW with convergence. No Γel/86eV Ω/BEC/CLASS claims.

## 2026-11-22 fixed-volume charge-constrained electron-pair heat

- [x] Derive exact finite-charge susceptibility C_e±|nQ /T³ = B−D²/S from Fermi electron/positron variances and Cauchy–Schwarz; verify positive heat and independently finite-different ρ(T±1e−4T) with μ recalculated to keep **physical net charge** constant. [Audit](../FIXED-CHARGE-PAIR-CAPACITY-AUDIT.md), [tests](../../tests/test_fixed_charge_pair_capacity.py).
- [x] **Negative:** at T/me=.025, illustrative η=6e−10, the electron-only fixed μ/T heat is ~1086 times the correct fixed-nQ electron heat, although photon heat still dominates total EM heat. Fast pair chemistry is a hypothetical closure, not established via annihilation rate.
- [ ] Verify new scientific CI and reproduce small/large y, quadrature/cutoff, independent FD; resolve thermal asymptotic T≪me for realistic physical baryon charge.
- [ ] **P0:** derive actual e−+e+↔γ+γ microscopic finite-T Pauli/Bose source and test forward/backward equilibrium, charge/energy/entropy, e± freeze-out, full photon+ν heat and UV-valid χe/eν before nonlinear χ+Bose FRW. No Γel/86eV Ω/BEC/CLASS claims.

## 2026-11-22 QED pair-annihilation microscopic partial gate

- [x] Independently derive vacuum QED e−e+→2γ spin-averaged CM differential amplitude, exact Dirac σv and lab electron-rest Heitler/Geant4 σ; verify β threshold and identical-final-photon factor (<8e−13 at five β values). [Report](../QED-PAIR-ANNIHILATION-CONTROL.md), [tests](../../tests/test_qed_pair_annihilation_control.py).
- [x] FD-initial collision **LOSS ONLY** R with correct Møller flux and spin degeneracy; compare Gauss radial/angular orders and low-T NR threshold. R/(n_e− n_e+ πα²/me²)=0.994852 at T/me=.025, .295344 at T/me=1. Not a physical net e± annihilation chemical rate.
- [x] **Negative:** exact full Bose/Pauli local quantum event factors cancel thermal e−+e+↔γ+γ at μe−+μe+−2μγ=0 (including nonzero opposite electric μ). Loss-only R cannot drive a correct equilibrium.
- [ ] Verify real CI for this stage and previous charge-constrained heat stage; add full 5D finite-T differential QED F−R with photon stimulation and inverse, independently check all number/charge/energy/entropy moments, λ/angle/temperature grids and finite-time chemical freezeout; add neutrino/UV-dark portal before full χe–BEST Bose FRW. No physical γel/BEC/86eV Ω/CLASS.

## 2026-11-22 QED full thermal 5D pair forward/inverse and heat ledgers

- [x] Construct actual 5D CM scattering photons Lorentz-boosted to lab with exact tree QED differential amplitude and identical γ factor, full fermionic Pauli e± and photon Bose forward+inverse factors. Exact thermal FD/BE common-T cancellation at electron chemical sum zero including opposite charge μ; positive entropy for chemical/temperature affinity. [Report](../QED-PAIR-FULL5D-AUDIT.md), [source](../qed_pair_full5d.py).
- [x] Independent microscopic normalization: empty γ initial state/absence of final enhancement reproduces separate exact vacuum σv×FD e± loss at T/me=.025,.05,.1 to ~1e−14; energies/charge and photon creation ledger consistent to rounding. [Tests](../../tests/test_qed_pair_full5d.py).
- [x] **Negative:** despite forward Bose enhancement +34.326% at Te/me=1, inverse exactly cancels thermal equilibrium; treating forward rate as NET equilibrium depletion violates quantum detailed balance. p/T cutoff scan nonzero; no claim of uniform numerical accuracy.
- [ ] Verify real CI for full 5D thermal stage, and prior vacuum QED stage. **Next P0:** nonthermal mode-resolved full strong source C_e±(p), Cγ(k), quantum Pauli/Bose, finite-volume positivity/energy/charge/entropy and stiff time evolution; include electron/proton charge neutrality, neutrino weak, thermal photon/QED corrections, gauge-valid χ1/eν portal before physical γel and 86eV Bose FRW freezeout. No ΩDM/BEC/CLASS.

## 2026-11-22 nonthermal e±γ full 5D QED weak collision

- [x] Generalize exact on-shell full-Pauli/Bose QED 5D forward-minus-reverse event flow to arbitrary positive bounded isotropic e−/e+ and γ occupation functions (no thermal FD/BE closure), with distinct charge/energy and positive entropy weak moments. [Audit](../QED-PAIR-NONTHERMAL-WEAK-AUDIT.md), [tests](../../tests/test_qed_pair_nonthermal_weak.py).
- [x] Independent thermal FD/BE reference +.02 affinity matches nonthermal routine ~5.89e−8 relative at Te/me=.1; genuine nonthermal electron packet ±.15 flips R, both directions entropy positive. Photon ±.2 bump can flip sign and entropy remains positive.
- [x] **Negative:** raw F−R may lose digits at Te/me=1 (~3e−6 relative vs stabilized), and photon bump nrad24→32 changes R ~1e−5. Uniform high-precision nonthermal QED operator not established.
- [ ] Check real CI; **P0:** robust log-ratio numerics, phase-space mode-resolved C_e−,C_e+,Cγ with exact microreversible deposition and Pauli/Bose bounded energy/charge-conserving time evolution; q/angle/cutoff/dt independent error budget, then ν/eγ SM, gauge-valid χ1 mediator, full χ+Bose FRW. No physical Γel/86eV Ω/BEC/CLASS yet.

## 2026-11-22 gauge-covariant lepton vector portal and neutrino heat

- [x] Identify explicit anomaly-free candidate U1_(Le−Lmu) with χ vectorlike singlet, lepton doublet+singlet currents, **mandatory active νe/νμ** for this vector e_L+e_R portal; independent gamma5 trace confirms one-helicity ν_L vector χ parity-symmetric amplitude equals massless Dirac-vector spin average. [Report](../LEPTON-PORTAL-NEUTRINO-HEAT-AUDIT.md).
- [x] Conditional full 5D Pauli contact EFT equal-T Qχ from two ν flavors vs electron e±: ratio 1086.37 at Te/me=.1, 1.185 at1, 1.00179 at10; finite cutoff/grid scans; no physical mediator/coupling assumed.
- [x] **Negative:** single SM bath Tν=Tγ can give even **wrong net χ heat sign** after photon reheating. Illustrative ideal split Tν/Tγ=(4/11)^(1/3) at Tγ=.1me, Tχ=.9Tγ leads Qν/Qe≈−512; this is NOT actual ν temperature history or 86eV DM prediction.
- [ ] Check CI and **P0:** finite Z' mass t-propagator and independent chiral/electron 5D normalization, UV/EW completion for actual χ1/χ2φ, realistic ν weak decoupling and e±γ entropy, nonthermal Pauli χ+ν+e kinetics and spectral γel, Bose χ→φ+χ2 and FRW/Ω/CLASS after full convergence. Experimental constraints essential.


## 2026-11-21 finite mediator EFT validity gate

- [x] Analytically derive exact spacelike t-channel vector propagator correction and both relative-error conventions; record U(1)' charge constraint for χ1–χ2–φ Yukawa in [note](../FINITE-MEDIATOR-CONTACT-VALIDITY-NOTE.md).
- [ ] **Immediate:** insert finite-mV propagator into independent full-Pauli χe/χν 5D heat; verify with separate t quadrature, angular/cutoff convergence, split baths, thermal-tail bound and actual charge/mediator constraints. Analytic bound is not numerical validation or physical gamma_el.

## Finite mediator angular diagnostic (21 Nov 2026)

- [x] Independently validate exact isolated propagator angular integrals, including transport weighting and extreme mass ratios (analytic vs 32768-step Simpson; local test). This is **not** a full physical χe/χν rate; [report](../FINITE-MEDIATOR-ANGULAR-CONTROL.md).
- [ ] Replace contact propagator in full 5D Pauli χe/χν collision, check independent t quadrature and thermal-tail convergence, gauge charge assignment, and real plasma heat history.

## Finite-mediator 5D first gate (22 Nov 2026)

- [x] Insert exact spacelike vector propagator factor into independently boosted full-Pauli 5D χ-lepton heat with full angular spin trace; retain old contact default.
- [ ] Verify real public CI for contact/equilibrium/finite-mV test; independent nonconstant-numerator angular quadrature and n/cutoff/mediator scans for massive e and chiral ν, including split baths. No physical gamma or 86-eV claim.

## 2026-11-22 finite mediator nonconstant spin trace angular gate

- [x] Independently derive exact t-integrated **nonconstant** vector spin trace with finite mediator, stable small-b series and adaptive quadrature tests; local isolated Jn validation ≤3.5e−14 relative. [Report](../FINITE-MEDIATOR-ANGULAR-TRACE-AUDIT.md). Prior full 5D finite mediator CI #44968592762 succeeded.
- [ ] Check CI for new trace tests; independently integrate full boosted final-state χe/χν Pauli heat over t/azimuth (cannot factor vacuum angular trace), n/cutoff/mV scans, split ν/e baths and thermal-tail bounds. Then physical finite SM plasma, gauge completion and nonlinear kinetic eigenmodes; no 86eV relic/BEC/CLASS.

## 2026-11-22 independent transfer-angle full-Pauli finite-mediator quadrature

- [x] Implement exact rational transfer coordinate in boosted 5D Pauli heat, retain full outgoing-state Fermi factors, vector trace and energy transfer; add comparison/equilibrium tests and report.
- [ ] Check new CI and independently scan mV/T, nrad/nangle/nphi and incoming/outgoing cutoff, massive e and massless chiral ν with split baths. If any discrepancy, retain numerical negative and repair before interpretation.
- [ ] Actual UV-valid χ1/χ2φ portal, finite SM chemistry/weak bath and nonthermal kinetic spectrum; no 86eV relic/BEC/CLASS claim.

- [x] Identify actual failed 5D transfer-angle CI [44979564709](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44979564709): Python 3.11/3.13 SyntaxError from literal `\\n` in source; repair in `99a15ea9`.
- [ ] Verify CI after syntax repair, resolve remaining failures; only then independent n/cutoff/mV thermal heat convergence, including negative cases.

## 2026-11-22 adaptive forward mediator angular gate

- [x] Previous repaired HEAD 25de5e55 passed real Scientific CI #44997020168.
- [x] Independent adaptive scipy polar reference and exact propagator angular identity; three local regression tests pass, negative mV=.03 48-node transfer error +1.05%; previous transfer audit prose factor-two corrected.
- [ ] Verify new CI. **P0:** integrate adaptive polar over full p,k,mu 5D physical e and active ν Pauli weak source with independent n/cutoff/mV/T and bath-splitting errors; then physical gauge/UV mediator+χ1χ2φ and nonthermal Bose FRW. No gamma_el/86eV Ω/BEC/CLASS claim.

## 2026-11-22 full 5D independent adaptive finite-mediator angular pilot

- [x] Upgrade fixed-event adaptive polar reference to full incoming Gauss p,k,mu 5D physical Pauli heat. At mV=.03 nrad=nmu=10, nphi16, angular Gauss32 transfer −9.559% and cosine −7.876% against adaptive reference on same incoming grid, both wrong far above 1%.
- [x] Repeat two incoming resolutions nrad=nmu=8 and10; inspect positive/negative quantum heat signs for mlep=0 and mlep=.5 plus exact common-temperature cancellation. Source, compact JSON and tests recorded. Global cutoff and time convergence NOT proved.
- [ ] Prior adaptive fixed-event CI #45006486258 passed; verify real GH CI for new full 5D stage. **P0:** massive e and chiral ν split baths; independently refine incoming radial, μ, azimuth, forward-polar n, energy cutoffs/mV; quantify absolute and signed heat uncertainties before physical e±γν+UV dark portal/Bose kinetic FRW. No gamma_el/86eV ΩDM/BEC/CLASS.

- [x] **Negative CI:** Scientific #45006858598 failed 1/812 due unresolved 12-node conventional Gaussian vs periodic azimuth (relative 1.569e−5). Independent nphi24 agrees ~1.39e−8. Updated regression retains explicit low-nphi alias and tests converged nphi24; new CI pending.
- [ ] Check corrected CI, continue independent full 5D nrad/nmu/nphi/cutoff convergence; preserve low-angle error report, then split ν/e and physical SM/UV + Bose kinetics.

## 2026-11-22 split lepton finite-vector heat with error cancellation

- [x] Conditional full Pauli 5D adaptive χ−e± and χ−two-active-ν with assumed Tν/Tγ=.713765855, Te=1,Tχ=.9,Mχ=.1,G=.001,mV=.03,.1,.3; Qe>0, Qν<0 under this assumption, total <0; separate heat errors can exceed 7-10% even if net relative error ~3%. [Report](../FINITE-MEDIATOR-SPLIT-LEPTON-AUDIT.md).
- [x] Negative radial scan at mV=.03: incoming nrad=nmu8→10 changes total heat by ~0.52%, so even polar-adaptive thermal total not fully converged; no physical cosmological sign claim.
- [ ] Check repair CI #45007251387 and new split-stage CI; P0 independent absolute channel q/polar/phi/cutoff convergence and physical e±γν weak history, viable gauge mediator and complete nonthermal Bose production. CLASS-ready false.

## 2026-11-22 real φ / χ1,χ2 gauge-consistency P0

- [x] Analytic U1' selection: qφ−q2+q1=0; original single real φ necessarily qφ=0, hence q1=q2. Nonzero diagonal χ1 vector portal enforces allowed χ2 vector coupling. A χ1-only charged-Z' and χ2-neutral renormalizable real-φ decay is forbidden in this minimalist completion. [Audit](../REAL-PHI-YUKAWA-GAUGE-CHARGE-AUDIT.md), 4 local rational tests.
- [x] Preserved alternative complex charged φ or mixed gauge-breaking completion as a NEW model, not reinterpretation of original real φ; Z2 remains consistent for equal-charge case.
- [x] Previous azimuth repair Scientific CI #45007251387 succeeded on 78662b.
- [ ] **Immediate P0:** check new gauge-claim CI and previous split-bath CI. Choose coherent χ1/χ2/φ U1' + Z2 theory. If real φ and equal gauge charges, derive missing χ2e/χν finite-T full Pauli collisions and simultaneous SM heat before physical Γel. For charged complex scalar, rederive all width/Bose degrees and stability. Then actual mediator constraints, SM neutrino decoupling and nonthermal Bose FRW. No 86eV Ω/BEC/CLASS claim.

## 2026-11-22 χ2 SM interaction forced by equal-charge real-φ model

- [x] Compute conditional independent full 5D thermal Pauli vector mediator heat separately for mother m1=.1me and daughter m2=.02me with real mφ=86eV threshold open. Same U1' charge and prescribed same thermal fχi: χ2/χ1 ν heat ratios 1.274 (incoming n6),1.277(n8),1.280(n10); massive e heat ratios ~1.08/1.07. Does NOT give daughter abundance or relic. [Audit](../CHI2-LEPTON-HEAT-GATE.md).
- [x] **Negative:** daughter lepton scattering is not suppressed by daughter mass alone in this benchmark, so χ1-only BGK or 5D heat cannot represent a full real-φ equal-charge model without joint χ1+χ2+SM source.
- [ ] Prior gauge-charge CI #45007616038 success; verify daughter-heat CI; P0 coherent χ1/χ2/φ gauge Lagrangian and bounds, full χ2/e/ν Pauli nonlinear collision with real fχ2(t,q), finite e±γν thermo & Bose Yukawa/FRW independent convergence before 86eV claims.

## 2026-11-22 real Zprime population and gauge-breaking gate

- [x] Exact mass thresholds for χ1/χ2 pair annihilation to 2 Zprime, Zprime decays to νe,νμ, to both dark χ pairs and electron pairs. Light mV=.03/.1/.3me have different open channels; vacuum widths depend on individual gprime and gchi rather than only product G mV².
- [x] Conditional massive-vector g=3 thermal Bose EOS independently checked against K1/K2, delta gstarR near 3 at T=me for mV=.03/.1/.3me. **Negative:** virtual t-exchange-only model can omit physical real vector population and heat; no population or Neff demonstrated.
- [x] Original real Z2-odd phi cannot act as charged Z2-even Higgs-breaking S; separate complex S or Stueckelberg alternative in coherent theory. Previous daughter heat scientific CI #45007896155 success.
- [ ] Check new CI. P0 true vector quantum Γ/H and population, separate χ1/χ2 Pauli scattering, neutrino/eγ/φ heat, gauge/BBN constraints and real Bose+FRW before 86eV Omega/BEC/CLASS.

## 2026-11-22 conditional Zprime-neutrino decay clock

- [x] Exact two-active-ν rest partial width, independent MB momentum integral K1/K2 thermal time dilation and radiation H at explicitly assumed gstar=10.75. With fixed arbitrary Gme²=.001 and conventional |gχqχ|<=sqrt4π, minimum gprime gives ΓMB/H=3.37e3,4.58e6,3.16e9 at mV=.03,.1,.3me. No lab viability or actual Zprime population follows.
- [x] Analytic ν FD μ=0 final-state Pauli bound: (1−fν)(1−fνbar)>=1/4. Only a conditional lower decay-loss clock for MB-distributed existing Zprime, NOT equilibrium or chemical relaxation.
- [ ] Check GH CI. **P0:** full moving Zprime↔νν finite-T mode-resolved Pauli/Bose inverse production with exact boosts and finite neutrino heat/entropy, then actual Γ/H and number-changing modes plus χ1/χ2 SM, Yukawa/Bose FRW. No physical Gammael/86eV Omega/BEC/CLASS.

## 2026-11-22 full moving Zprime-neutrino quantum 1<->2 inverse production

- [x] Derive exact finite-T Zprime ννbar 1<->2 quantum loss/inverse with full Pauli final ν, BE vector occupation and on-shell moving-vector angular E1,2=(E±p cosθ)/2; show A(p,T) closed log and p→0 tanh(mV/(4T)), independent Gauss test, equal-T BE fixed point, entropy and energy/number ledgers. Not a χ elastic operator.
- [x] In explicitly fixed thermal ν reservoir T=me, zero ν chemical potentials, initial empty mediator and numerical product Gme²=.001 with |gχeff|=sqrt4π, compute physical restricted inverse vector fill/H at mV=.03,.1,.3me =1.34973e3,1.84904e6,1.33809e9; independent radial convergence. This is conditional reservoir-only, not cosmological thermalization.
- [x] **Negative:** MB unblocked clock differs by ~2.4× from fully quantum FD/BE inverse filling; neither arbitrary elastic BGK gamma nor prescribed equilibrium vector gstar justified for real FRW without backreaction.
- [ ] Verify GH Scientific CI for new Zprime quantum stage and previous MB clock. **P0:** coupled finite ν energy/entropy and vector b(p) time/H, additional vector/dark/e± channels, then both χ1/χ2-SM Pauli strong plus original φ Bose/Yukawa/FRW, full independent q/dt/IR/observational gates.

## 2026-11-22 closed Zprime/ν number+energy charge equilibrium counterexample

- [x] Exact reduced charge for ννbar↔Zprime only: Nν+Nνbar+2NV and total E; with FD ν/antiν common number mu and BE vector muV=2mu solve fixed-volume maximum-entropy equilibrium independently with n128/256 radial and formal massless-3-pol fugacity sums. Initial bV=0, νmu=0 leads final muν/T≈−1.136 and nV only ~0.17 of zero-mu Bose T0, entropy +7.3% at mV/T0=.03..3.
- [x] **Negative:** a large quantum inverse fill/H does NOT prove final ν and Zprime chemical potentials vanish. In reduced system an implicit zero-mu thermal reset breaks conserved Nν+2NV and energy; full weak ννbar↔e+e− can invalidate reduced conservation and must be modeled.
- [ ] Check GH CI. **P0:** positive finite ν reservoir + exact vector mode time integration with T,mu determined by N/E, entropy, q/dt convergence; then real ν-e±γ weak number changing and vector/χ1/χ2+φ Bose/FRW. No physical Gamma_el, 86eV Omega/BEC/CLASS.

## 2026-11-22 finite ν chemical/thermal reservoir and Zprime Bose 1<->2 time

- [x] Evolve Zprime b(p) via exact mode exponential of physical ν Pauli+Bose inverse source; update finite massless ν FD(T,μ) to preserve true reduced Nν+2NV and total E, recompute entropy each step. Under explicit instant ν self-rethermalization closure, τ5000 steps500 reaches independent constrained equilibrium T to 7e−8 and µ/T to 1e−7, N/E errors ~1e−14, entropy monotonic. [Audit](../VECTOR-NEUTRINO-FINITE-BATH-TIME-AUDIT.md).
- [x] **Negative numerical:** τ100/20,40,80 has visible ~first-order splitting drift despite exact charge/energy/entropy; τ3000 in 20 substeps overdraws ν number budget and is rejected. Convergence and stiff adaptation remain P0.
- [x] Prior quantum vector finite-T CI #45009062951 on 9582b5a success.
- [ ] Verify new CI and closed-equilibrium CI; **P0:** genuine nonthermal ν strong source and SM weak ννbar↔e± number-changing, finite e±γν entropy and Zprime dark pair channels; gauge constraints and both χ1/χ2 full Pauli + real-φ Bose FRW, no physical Γel/86eV ΩDM/BEC/CLASS claim.

## 2026-11-22 nonthermal ν strong mode quantum source

- [x] Derive independent ν-mode Cν(k) from on-shell moving unpolarized vector decay with fν=fνbar arbitrary isotropic, exact Pauli/Bose and 2gV/gν normalization; independently integrate vector weak p,cosθ and ν strong k,p and compare number, energy, entropy. Four nonthermal shapes pass to <1e−9 at mV=.3 n64, six local tests pass. [Audit](../VECTOR-NEUTRINO-NONTHERMAL-STRONG-AUDIT.md).
- [x] **Negative:** exact inverse support pmin(k)=abs(k−mV²/4k) has cusp at k=mV/2. Unsplit 32-point k quadrature causes artificial ~−3.1e−4 weak/strong number discrepancy, reduced by splitting; mV=.03 needs far finer grids; same-cut agreement does not establish pmax tails.
- [x] Previous closed ν/V finite bath time CI #45009845364 success.
- [ ] Verify new CI. **P0:** coupled fν(k,t), fνbar(k,t), bV(p,t) physical positive strong time integration with actual Pauli/Bose gain/loss, exact N/E/entropy and q/dt/cutoff, no instantaneous ν-FD reservoir. Then weak e±γ/ν physics, gauge viable Z′/χ1/χ2/φ Bose and full FRW observables; Γel/86eV Ω/BEC/CLASS unproven.

## 2026-11-22 distinct nonthermal ν/antiν/V bounded event time solver

- [x] Independent full quantum 1↔2 positive event graph on exact integer energy shells with separate ν and antiν FD occupations, vector BE, implicit bounded extent and eventwise exact reduced number, lepton asymmetry and energy, entropy monotonic; local numeric pilot and permanent stage/JSON.
- [x] **Negative:** single energy lattice angular quadrature can be badly biased: weak source N −21% at h=.4 vs 144×144 continuous reference, h=.05 still +0.38%; Emax4→8 large tail; dt1→.125 sizable low-energy bV shift. Do not call it continuum strong time solver or physical dark matter.
- [x] Prior syntax repair 64c394 has real Scientific CI #45010866936 success.
- [ ] Verify new GH CI. P0 fix angle/energy quadrature without sacrificing microreversibility, construct convergent higher-order time and multi-shell strong source; then finite ν/eγ and gauge-valid χ1χ2/φ Bose FRW. CLASS-ready false.

## 2026-11-22 exact angular shell overlap diagnostic

- [x] Shift ν and antiν energies to half-integer cells while maintaining EVector=jh=ki+kl. Use true length of each daughter's energy cell within k∈[(E−p)/2,(E+p)/2] for microreversible event rate W. Flat CM angular probability sums to unity per mother shell to <4e−15 at four grids; no posthoc rate renormalization.
- [x] **Positive scoped:** at h=.2 same daughter grid center-hit vs overlap, number weak rate relative to independent continuum improves +5.2037%→+0.0824%, entropy +5.3046%→+0.8059%. [Audit](../VECTOR-NEUTRINO-SHELL-CELL-OVERLAP-AUDIT.md).
- [x] **Negative:** exact angular weights still leave 2.229% energy-source bias at h=.2, 0.592% at h=.05, and unresolved Emax/dt errors. No strong continuum quantum cosmology rate claim.
- [ ] Check both real CI runs; P0 vector energy/radial projection and higher-order positivity/entropy time, then finite ν/eγ/χ1/χ2/φ gauge/FRW and only later 86eV observables.

## 2026-11-22 vector E-cutoff endpoint density measure gate

- [x] Derive exact positive vector shell number density integral from E p(E)dE, clip upper Emax and lower mV; compare full h boundary, upper trapezoid and exact cell on identical energy/angular quantum event graph with 3+ h. [Audit](../VECTOR-NEUTRINO-SHELL-RADIAL-ENDPOINTS-AUDIT.md).
- [x] Positive scoped: h=.05 energy weak bias +0.592%→−0.032%, angular normalization and discrete eventwise invariants unchanged. Prior reversible-event full Scientific CI #45011435812 success.
- [x] **Negative:** h=.05 number −0.142%, h=.4 number error worsens −0.712→−2.694%; missing E=mV threshold cell, time splitting and cutoff not solved. No physical Γelastic/86eV Omega/BEC/CLASS.
- [ ] Verify CI for angular stage 174e9b (#45011770137 pending last check) and new radial stage; P0 physical exact cell energy and threshold, high-order dt and continuum sources before SM/χ/Bose FRW.
- [x] Real Scientific verification of angular overlap 174e9b: [run #45011770137](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45011770137) succeeded; vector endpoint/new stage CI still to confirm.

## 2026-11-22 exact reaction-network conserved-charge topology gate

- [x] Analyze stoichiometric integer event incidence for Zprime↔nu+antinu with two independent modular primes and SVD; C/L/E left residual zero, observed original integer grid nullity5 (2 artificial) vs shifted angular-overlap nullity3 (physical) for h=.4,.2,.1 at mV=.4,Emax4. [Audit](../VECTOR-NEUTRINO-EVENT-INVARIANTS-AUDIT.md).
- [x] **Negative structural, not mere quadrature:** original integer grid has two isolated top k=Emax ν/antiν bins, whose occupation is exactly frozen by every event for any timestep. H-theorem/CLE conservation cannot prove physical constrained equilibrium on that grid.
- [ ] Check new GH CI and prior stage. **P0:** independent maximum-entropy FD/BE constrained by actual discrete C,L,E with ν−antiν asymmetry; long-time event dynamics, full nullspace/slow modes and h/dt/Emax, then physical weak ν–SM, χ1χ2φ Bose FRW. CLASS-ready false.
- [x] Previous shell-radial endpoints commit b37b7dd passed real [Scientific CI #45012189092](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45012189092); new stoichiometric topology stage CI pending.

## 2026-11-22 independent maximum entropy versus nonthermal strong time

- [x] Independent three-constraint maximum-entropy FD ν, antiν / BE vector with nonzero lepton asymmetry and true shifted discrete C,L,E, at h=.4,.2,.1. Compare to actual time evolution preserving all ν,antiν,bV distributions without thermal reset. τ500(dt1) h=.2 mode abs difference .001037, entropy gap1.86e−6; τ1000 .000339, gap3.99e−7.
- [x] **Negative:** τ100 h=.2 vector lowest bV for dt2/1/.5/.25=.75305,.743975,.746678,.748743, not monotone, and h=.4→.1 equilibrium T/μ changes. No continuum physical time-rate certificate.
- [x] Previous nullspace proof baa501e full Scientific CI #45012556392 successful.
- [ ] Verify new CI. P0 compare symmetric / reverse event order vs actual higher-order time and fully implicit, with 3+ dt/h/cutoff spectral checks, then ν-e±γ/χ1χ2φ Bose FRW. No Γel/86eV relic/CLASS claim.

## 2026-11-22 exact Pauli/Bose reaction flow and Strang temporal audit

- [x] Derive exact local quantum 1↔2 Riccati reaction extent and compare 18 independent event/time ODE integrations at abs≤3.2e−16; use palindromic forward/reverse half-events, no thermal reset, hard quantum bounds, conserve discrete C/L/E and grow entropy.
- [x] Positive scoped: h=.2 full three-species nonthermal spectral time error dt=.0625/.03125/.015625 vs dtref=.00390625 = .007275/.002098/.000522, ratios 3.47/4.02, consistent with formal order2 on fixed graph.
- [x] **Negative:** h=.1 is stiff at tested dt and does not show order2, errors .02321/.01417/.00722; previous angular/radial/number/cutoff issues remain. This is **NOT** physical χ Gammael or 86eV relic production.
- [ ] Verify current GH CI; P0 h=.1 smaller dt and 3+ Emax/h/shape comparison, then true finite e±γν/χ1χ2φ SM/FRW Bose, no CLASS-ready.

## 2026-11-22 quantum Strang strong-time stiffness and self-convergence

- [x] Independent full-three-spectrum dt vs dt/2 self-convergence for mV=.4 Emax4, tau1, h=.4/.2/.1, dt=1/8…1/1024. On coarse h=.1 first relative error ratio .794 despite decreasing dt, proving non-asymptotic strong-event stiffness. [Stage](../../archive/stages/2026/11/2026-11-22-vector-neutrino-strang-dt.md).
- [x] Positive scoped: refined h=.1 ratios approach 3.55,3.876,3.968 while C/L/E ~1e−15, H entropy and Pauli/Bose. This validates second temporal order only on finite quantum reaction graph, not physical continuum or χ–SM.
- [ ] Verify parent/new GH CI. P0 adaptive full-spectrum error-step controller, source h/Emax/shape independent convergence and gauge/SM χ1/χ2φ Bose+FRW. No Γelastic/86eV relic/CLASS.
