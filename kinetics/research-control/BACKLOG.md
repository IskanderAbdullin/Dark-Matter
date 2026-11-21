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
