# Stage 2026-11-22 — exact log-transfer outgoing CM angle for finite-vector χℓ Pauli 2→2

- **ID:** DM-KINETICS-20261122-CHI-LEPTON-LOG-TRANSFER-CM
- **Session UTC start:** 2026-11-22T14:52:42Z (independent clock). Exact substage end not separately recorded.
- **Preceding research HEAD:** [78d401cdba2956c61aa664e1f30b44e2e9c9447e](https://github.com/IskanderAbdullin/Dark-Matter/commit/78d401cdba2956c61aa664e1f30b44e2e9c9447e), [actual Scientific verification #45036745455 success](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45036745455). Draft PR #1 unchanged.
- **Project context:** [GOAL](../../../../GOAL.md) → [PLAN](../../../../PLAN.md) → [kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md) → [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md) → [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Research hypothesis, negative alternatives, criteria

Momentum segmented source reduces one radial alias, but CM cosθ Gauss still nonmonotone, particularly when mV² ≪2p*². Fix physical mediator propagator² angular integration **without introducing a phenomenological relaxation rate**. Validate against independent exact angular integrals I0 and transfer-weighted I1 for A=2p*²/mV² from .01 to 10000; test a known defective I0-only compensation, and two full weak 5D Pauli synthetic regimes with different mediator masses and 3+ angle resolutions. Do not discard untransformed negative outcomes or assume absolute physical convergence.

## Derivation and physical meaning

t=−2p*²(1−c). Define z=1−c and A=2p*²/mV². Propagator²=(1+Az)^−2. Analytic tests: I0=2/(1+2A); I1=[log1p(2A)−2A/(1+2A)]/A². Transform y=log1p(Az), y_max=log1p(2A), c=1−expm1(y)/A; exact phase-space product propagator² dc=exp(−y)/A dy. Gauss y gives nonnegative weights; original on-shell χ+ℓ→χ′+ℓ′ boost, matrix-element NUMERATOR, Pauli F/R and microscopic entropy are unchanged. Alternate variable u=(1+Az)^−1 yields exact I0 weight but **does not** stabilize I1 automatically.

[Independent analytic+full 5D report](../../../../kinetics/CHI-LEPTON-LOG-TRANSFER-ANGLE-AUDIT.md), [source](../../../../kinetics/chi_lepton_outgoing_log_angle.py), [modified independent full event with default uniform mode preserved](../../../../kinetics/independent_chi_lepton_5d_quadrature.py), [tests](../../../../tests/test_chi_lepton_outgoing_log_angle.py), [data](../../../../kinetics/results/chi-lepton-outgoing-log-angle-2026-11-22.json). Reproduce using OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_chi_lepton_outgoing_log_angle.py; python -m kinetics.chi_lepton_outgoing_log_angle --full. Locally **10 tests passed**; new CI status not yet independently confirmed.

## Verified positive and negative numerical evidence

At A=1000,ncos8 uniform c I0 error −93.1%, I1 −19.8%, u-cancel I1 −33.9%; log-transfer I0 ≈−4.4e−10 and I1 ≈+6.6e−11. At A=10000 even the u-only I1 misses by −50.3%, while log-transfer I0 error −1.42e−8 and I1 +1.60e−9. This is a strict independent mathematical kernel check.

On full five-dimensional χℓ Pauli weak Q with synthetic smooth Fermi and **narrow mV=.01** (Mχ=.1,mℓ=0,segmented radial8×4,nμ6,nφ12,cut12), outgoing ncos8/16/32/64 log-transfer Q=(1.7999415280,1.7999394500,1.7999394481,1.7999394481)e−17 while uniform ncos8/16/32/64 Q=(.88454130379,1.140027874,1.387397025,1.600559905)e−17. Uniform ncos32 still undercounts mapped ncos64 by **22.92%**; mapped ncos8 differs by only ~0.000116% at same finite incoming measure. Pauli entropy behaves similarly.

For oscillatory f and **mV=.1** (segmented radial12×5,nμ8,nφ16,cut12), mapped ncos8/16/32/64 Q=(4.674999,4.662545,4.658706,4.658669)e−15 and H=(7.432199,7.432370,7.432384,7.432386)e−13. Mapped 32→64 Q shift ~0.00081%. But uniform cosθ Q over coarse angular grids is strongly nonmonotone. These values are scoped **synthetic conditional 5D** improvements, not full physical rates.

## Negative scientific boundary and next step

Exact I0/I1 and dramatic forward-mapping improvement do NOT establish convergence in incoming radial, relative μ, outgoing φ, UV/cutoff or lepton mass, nor viability of the gauge model. No full strong Cχ1,Cχ2,Cℓ with actual fχ2 from χ1→χ2+real φ Bose, no e±γν/Zprime thermal plasma, no FRW H or freezeout. No 86eV ΩDM, BEC, or CLASS-ready claim. Next executable stage: joint grid/cutoff/finite-mass error controls **with this same mapped physical amplitude**, chemical detailed balance and explicit true-gspin/SM channel bookkeeping. Preserve earlier negative heat and angular results in evidence register.

## Actual CI (append-only)

- [Scientific verification #45037526409](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45037526409) on 1e5b497a04c652a20b43cbdddb77efc697aa035b completed success. It certifies scoped mathematical/numeric CM angular controls only.
