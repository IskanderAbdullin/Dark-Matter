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
