# Stage 2026-11-22 — χℓ 5D segmented energy-guided momentum integration

- **ID:** DM-KINETICS-20261122-CHI-LEPTON-SEGMENTED-RADIAL
- **UTC:** active research run began 2026-11-22T14:52:42Z; exact per-stage elapsed time not measured, do not infer from task counts.
- **Previous branch SHA:** [c6d87db35a00bbce43af9be0885b4035b350139a](https://github.com/IskanderAbdullin/Dark-Matter/commit/c6d87db35a00bbce43af9be0885b4035b350139a), [real Scientific CI #45035908197 success](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45035908197). Draft PR #1 remains in public research/kinetic-equilibrium-audit-20261121; new stage CI pending.
- **Context:** [GOAL](../../../../GOAL.md), [PLAN](../../../../PLAN.md), [specific kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md), [constraints](../../../../kinetics/RESEARCH_GUARDRAILS.md), [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Hypothesis and numerical falsification

The previous physically allowed oscillatory χℓ Pauli weak f makes incoming unsegmented p,k Gaussian 5D heat flip sign at small nrad and cutoff. Test energy-guided *momentum* Gauss slices versus energy-coordinate Gauss near the massive square-root threshold; separately vary nlocal and nsegments before interpreting cutoff. Require physical radial number measure ∫p²dp exact to roundoff, finite microscopic entropy, explicit Q and H, and retain any angular residual rather than asserting global continuum convergence.

For each species define E_j=m+j ΔE and p_j=√(E_j²−m²), with nlocal Gauss points in [p_j,p_{j+1}], ΔE≤1 (arbitrary Tref=1). The p-space cell volume is polynomial, sum of weighted p² is pmax³/3 to 3e−14, unlike direct E-Gauss threshold cusp (m=.5, nseg12,nlocal5 relative volume defect +4.52e−7). Exact scattering, full F/R Pauli and Qχ=-Qℓ follow the independent [5D source](../../../../kinetics/independent_chi_lepton_5d_quadrature.py); no g/lepton flavor physical calibration or χ2 dynamic yet.

## Reproducible inputs and data

[Code](../../../../kinetics/chi_lepton_segmented_5d_gate.py), [modified reference with default uniform-p preserved](../../../../kinetics/independent_chi_lepton_5d_quadrature.py), [tests](../../../../tests/test_chi_lepton_segmented_5d_gate.py), [full audit](../../../../kinetics/CHI-LEPTON-SEGMENTED-5D-AUDIT.md), [JSON data](../../../../kinetics/results/chi-lepton-segmented-5d-2026-11-22.json).

Reproduce locally:
    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_chi_lepton_segmented_5d_gate.py tests/test_chi_lepton_quadrature_gate.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.chi_lepton_segmented_5d_gate --full

**Result:** eight focused tests passed locally in 24.3s; no RNG/seed. Synthetic Mχ=.1,mℓ=0,mV=.1,G=.001, historical numeric gχgℓ=4×4, nμ8,ncos32,nφ16, two nonthermal prescribed f forms. This is NOT physical Γel/H.

## Actual results and negative uncertainty

At cutoff12 and ΔE≤1, nsegments12*nlocal3/4/5/6 yields oscillatory Qχ=(4.645350,4.499933,4.568756,4.556594)×10−15 and positive H≈7.53–7.54×10−13. Compare unsegmented uniform-p nr24 Qχ=−1.398868×10−14 (wrong sign) and unsegmented nr64 Qχ=+4.556191×10−15. Segmented12×6 versus uniform64 agree to ~0.009% **at the same fixed angular/cutoff setting**. However segmented8×3 Qχ=2.27825×10−14 and segmented8×4 Qχ=1.29928×10−15, so segmentation alone does not guarantee convergence or high-order resolution.

When raising cutoff **with ΔE kept ≤1** (and nlocal5), oscillatory Qχ for cutoff8/12/16/20 = (5.365358,4.568756,4.594670,4.594453)×10−15. Cutoff16→20 shifts Q ~0.0047%, H ~0.0034%; this is a deliberately scoped cutoff check. Smooth f Qχ=(8.75901,9.06987,9.08250,9.08291)×10−14. The cutoff8 estimate differs materially, so all three-or-more cutoff data retained.

**Still negative:** for same segmented12×5 with cutoff12, outgoing CM polar orders16/24/32/48/64 produce Qχ=(4.87480,4.60634,4.56876,4.62164,4.65110)×10−15 and entropy (6.8560,7.3998,7.5405,7.5429,7.4984)×10−13. Q and S remain nonmonotone under angular refinement. A physical sub-percent Γelastic **cannot** be inferred from the stable radial volume/UV alone; fully strong χ and lepton mode operator and real φ/SM/FRW remain absent.

## Decision and executable follow-up

**Partial progress:** explicit energy-guided momentum segments and fixed-resolution cutoff comparisons resolve one demonstrated source of incoming radial numerical alias; historic unsegmented sign-flip remains a true negative. **Not established:** continuum angular/radial combined Q/S uncertainties, field-theoretic portal rate, physical χ1/χ2+real φ and SM bath. Next P0: forward-peaked CM angular segmentation or adaptive transform, cross-check independent 3+ ncos/nμ/nφ/radial/cutoff; then physical full Pauli+gauge+FRW. No 86eV Ω/BEC/CLASS claims. Update STATE/EVIDENCE/BACKLOG/NEXT_ACTION, commit+push, proceed if time budget allows.
