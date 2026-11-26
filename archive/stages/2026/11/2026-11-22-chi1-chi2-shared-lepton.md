# Stage 2026-11-22 — first nonlinear χ1+χ2 with one shared Pauli lepton reservoir

- **ID:** DM-KINETICS-20261122-CHI1-CHI2-SHARED-LEPTON-STRONG
- **Session start UTC:** 2026-11-22T14:52:42Z; test independently rerun at 15:30:08Z; no invented exact stage duration.
- **Previous HEAD/PR:** [b08c4cae4da2d8d6cac362c13719cceee8791bb4](https://github.com/IskanderAbdullin/Dark-Matter/commit/b08c4cae4da2d8d6cac362c13719cceee8791bb4), draft PR #1, real [Scientific CI #45038184612 success](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45038184612). This stage's GH CI is to be checked.
- **Context:** [GOAL](../../../../GOAL.md) → [PLAN](../../../../PLAN.md) → [kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md) → [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md) → [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Research hypothesis and independently testable criterion

Previously χ1 and χ2 scattering were numerically evaluated as two separate full Pauli 5D WEAK channels, not as **two evolving χ occupations on a single evolving lepton phase-space distribution**. Build physical elastic χaℓ↔χaℓ event source from independent on-shell 5D spin trace (full F/R Pauli) with **the same logit reconstruction/deposition weights** and correct single lepton grid shared by both channels; verify χ1 and χ2 individual elastic number and lepton number conserved, Qχ1+Qχ2+Qℓ=0, quantum entropy nonnegative, three independent chemical μ at common T fixed point, and adaptive Fermi-bounded nonlinear time, without any artificial FD resets/clipping.

## Formulas, source, and independent checks

Synthetic Mχ1=.1,Mχ2=.02,mℓ=0 (and fixed point mℓ=.5),mV=.1,G=.001,Tref=1, **coarse** incoming 5D nr4,nμ4, outgoing mapped CM8×φ8, cut12, gspin4 each inherited (not physical ν/e flavors). One common lepton grid covering both χ energy ranges, two separate χ grids with correct mass thresholds, exact radial zeroth cell weights w=g/(6π²)[(U²−m²)^(3/2)−(L²−m²)^(3/2)]. Microscopic 2→2 elastic event uses F=fχifℓ(1−fχi′)(1−fℓ′), R=inverse and deposits each signed χ channel into its own nodes but *both lepton updates into the SAME lepton nodes*. Matching barycentric logit(Fermi) reconstruction gives Sdot_node=∑W(F−R)ln(F/R)≥0, each species N and E invariant for pure elastic.

[Full strong source and time pilot](../../../../kinetics/chi1_chi2_shared_lepton_pauli.py), [four tests](../../../../tests/test_chi1_chi2_shared_lepton_pauli.py), [full audit](../../../../kinetics/CHI1-CHI2-SHARED-LEPTON-PAULI-AUDIT.md), [machine results](../../../../kinetics/results/chi1-chi2-shared-lepton-pauli-2026-11-22.json). Existing [independent weak](../../../../kinetics/independent_chi_lepton_5d_quadrature.py) extended to accept a PAIR OF CALLABLE Fermi functions, retaining historical synthetic labels. Reproduce: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_chi1_chi2_shared_lepton_pauli.py. Four local tests passed at 15:30:08Z.

## Results, historical negatives, and precise claim boundary

Initial current-distribution source: Qχ1=+8.469023562689831e−14, Qχ2=+1.0384403681255402e−13, ONE common lepton Qℓ=−1.8853427243945206e−13; entropy source direct3.4835786974671425e−14 vs same nodal-logit source3.483578697467172e−14. Nχ1/Nχ2/Nℓ source residuals ~1e−29 to1e−28. Independent earlier 5D WEAK evaluator for each Mχ and same fℓ matches BOTH χ heat and combined entropy to <2e−11 relative. Common-T μχ1=.1, μχ2=−.05, μℓ=−.08 with mℓ=.5 yields microscopic stationary source.

Real nodewise time simulation (not cosmological) to numerical τ=1e10 with adaptive first-order two-half Euler, atol1e−7: 11 accepted,7 rejected, smallest positive step entropy+2.1139e−05,total entropy+0.000340048; Nχ1/Nχ2/Nℓ drifts <5.6e−17, combined E drift <4.5e−16, all Fermi f strictly interior (min χ1 2.35e−12, χ2 8.26e−14, lepton3.42e−11). Finite Fermi cell capacities sum exact analytic radial density within 2e−14.

**Nontrivial negative/limitations:** both χ species are populated by independent **synthetic** initial functions; the χ2 spectrum is not physically produced by χ1↔χ2+φ. Historical 4×4 degeneracy does not represent physical SM neutrino/e± spin/flavor. The 5D nrad4/CM8 weak rates are **already disproved** as quantitatively reliable by earlier up-to-55× oscillatory Fermi convergence negative. Adaptive local tolerance does NOT certify global time/freezeout error and is first order. The individual χ numbers will cease to be invariants when Yukawa conversion opens; currently elastic-only. No real φ Bose, Zprime number, SM e±γν weak chemistry, H/FRW, physical Γ_el, ΩDM, BEC or CLASS-ready.

## Next executable scientific step

Verify new GitHub Scientific CI and preserve failures. Priority: derive real neutral φ Yukawa χ1↔χ2+φ fully quantum Bose/Pauli on-shell decay/inverse source, threshold Mχ1>Mχ2+mφ, exact chemical equilibrium and total fermion number/energy, partner-mode occupation; then actually couple to the tested shared finite lepton strong system and independently converge physical 5D q/angle/cutoff at ≥3 grids, correct flavor/couplings and full SM/FRW. No terminal scientific result claimed.

## Scientific CI status (append-only)

- Previous χ1+χ2 with a single shared lepton nonlinear Pauli stage [run #45038900605](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45038900605) status in_progress / conclusion null at this checkpoint; no physical φ produced by that stage.
