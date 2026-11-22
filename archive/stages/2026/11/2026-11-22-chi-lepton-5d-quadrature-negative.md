# Stage 2026-11-22 — independent χℓ 5D radial/angular/cutoff weak source nonconvergence

- **ID:** DM-KINETICS-20261122-CHI-LEPTON-QUADRATURE-NEGATIVE
- **Date:** 2026-11-22 UTC. Active research session started **2026-11-22T14:52:42Z** by direct container clock; no per-stage exact end time asserted.
- **Branch / prior SHA:** research/kinetic-equilibrium-audit-20261121, [12e3874c](https://github.com/IskanderAbdullin/Dark-Matter/commit/12e3874ccc50a89609add0a7b6781b2d71f73d88), draft PR #1. Its Scientific CI [#45017550678](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45017550678) independently confirmed success. Stage's own CI pending.
- **Goal routing:** [GOAL](../../../../GOAL.md) → [PLAN](../../../../PLAN.md) → [kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md) → [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md) → [state](../../../../kinetics/research-control/STATE.json)/[next](../../../../kinetics/research-control/NEXT_ACTION.json).

## Hypothesis and falsifying test

Former χℓ finite-vector exact local energy, Pauli and entropy gates do not guarantee that a low-resolution 5D weak heat Qχ is physically reliable. Repeat the same real on-shell amplitude and synthetic nonthermal Fermi spectra in an independent NumPy/SciPy implementation, compare against previous FD-specialization baseline, then scan incoming radial, incoming angle, outgoing CM polar/azimuth and kinetic cutoff **one at a time** and several joint high resolutions. Success means actual numerical stability on independent axes, not merely positive entropy.

## Method and reproducible evidence

Independently coded p/k/mu × CM polar/azimuth vector scattering quadrature, relativistic Mχ=.1, mℓ=0, mV=.1, G=.001, historical gχgℓ=4×4, synthetic Tref=1. Exact per-event E+ω=E′+ω′; F=fχfℓ(1−fχ′)(1−fℓ′), R=inverse and H=(F−R)ln(F/R). See [full report](../../../../kinetics/CHI-LEPTON-5D-QUADRATURE-NEGATIVE-AUDIT.md).

- [Independent evaluator](../../../../kinetics/independent_chi_lepton_5d_quadrature.py).
- [Five-axis scanner](../../../../kinetics/chi_lepton_5d_quadrature_gate.py).
- [Four tests](../../../../tests/test_chi_lepton_quadrature_gate.py).
- [Compact results](../../../../kinetics/results/chi-lepton-5d-quadrature-gate-2026-11-22.json).
- Locally executed: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_chi_lepton_quadrature_gate.py → **4 passed**.
- Full reproduction: OPENBLAS_NUM_THREADS=1 python -m kinetics.chi_lepton_5d_quadrature_gate --full.

## Actual positive and negative results

Independent synthetic FD Tχ=.9,Tℓ=1 baseline matches old weak Q=6.935904667187969×10^-14, S=7.706560741319965×10^-15; at common T with separately nonzero chemical potentials μχ=.2, μℓ=−.1, Q=2.87×10^-28, S=2.02×10^-42. Both nonthermal shapes have microscopic per-event entropy ≥0, energy defect <5×10^-14, Qχ+Qℓ=0.

**Negative:** At fixed nμ8,ncos24,nφ16,cutoff12, oscillatory Qχ with nrad=4,6,8,10,12,16 is (30.9727,−6.4605,+2.0409,−41.1258,+14.2008,+0.0091)×10^-14, with multiple **sign flips**. At nrad8,nμ8,ncos24,nφ16 cutoff=8,12,16,20 yields (−35.1256,+2.0409,−25.1233,−27.6789)×10^-14. For smooth f, outgoing cosθ=8→32 moves Q from 7.6719 to 9.0674 ×10^-14 (~18.2%).

Separate multi-axis refinement for oscillatory f: former inherited nrad4,nμ4,ncos12,nφ12 gives Q=2.55856732856×10^-13, versus nrad64,nμ16,ncos64,nφ32 Q=4.614781070×10^-15 (**55.44-fold smaller**). Refined conditional estimate is not a continuum certificate: ncos24..96 at fixed nrad40,nμ12,nφ24 still changes Q and H. All figures are synthetic numerical tests; no physical Γel or freeze-out.

## Scientific decision and next gate

Old 10-bin logit Pauli time pilot passed consistency but its underlying low-resolution 5D rates are not quantitatively physical. **Do not overwrite/discard those results:** interpret as structural/instrument only. P0: segmented/transformed incoming momentum integration to resolve nonthermal oscillations, adaptive outgoing forward angular quadrature, 3+ independent resolutions/cutoffs and Q/S error budget; only then physical χ1+χ2 real φ Bose+ν/eγ/Zprime/FRW. No χ Γel, 86eV ΩDM, BEC or CLASS-ready claim. Not terminal.
