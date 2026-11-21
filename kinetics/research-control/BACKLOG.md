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
