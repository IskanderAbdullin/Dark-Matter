# Stage 2026-11-22 — conditional two χ masses and finite lepton mass in Pauli 5D source

- **ID:** DM-KINETICS-20261122-JOINT-CHI1-CHI2-LEPTON-MASS-WEAK
- **Active session started:** 2026-11-22T14:52:42Z actual UTC clock. Exact per-stage duration not measured.
- **Starting HEAD/branch:** [1e5b497a04c652a20b43cbdddb77efc697aa035b](https://github.com/IskanderAbdullin/Dark-Matter/commit/1e5b497a04c652a20b43cbdddb77efc697aa035b), research/kinetic-equilibrium-audit-20261121, draft PR #1. New stage CI to verify.
- **Scientific context:** [GOAL](../../../../GOAL.md), [PLAN](../../../../PLAN.md), [special kinetic spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md), [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md), [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Predefined physical question and numerical criteria

In a neutral real-φ Yukawa χ1↔χ2+φ vectorlike equal-charge branch, BOTH χ species must scatter with SM leptons. Existing on-shell 5D finite-vector Pauli tests were mainly for Mχ1=.1 and massless ℓ; they do not show that the χ2 channel or charged-lepton mass works or is numerically stable. Check Mχ=.1/.02 and mℓ=0/.5 and two distinct admissible synthetic Fermi shapes with exact energy/entropy, independent FD chemical fixed point, and two joint quadrature levels for all eight cases; third joint level for selected channels. This experiment does not produce χ2 from χ1 decay or establish physical lepton flavor populations.

## Reproduction and explicit quantum model

Full relativistic on-shell finite-mV χℓ↔χℓ spin trace with forward/reverse Pauli, positive H=(F−R)ln(F/R), energy E+ω=E′+ω′; energy-guided momentum slices and log-transfer CM variable (previous exact I0/I1). Parameters Mχ=.1/.02,mℓ=0/.5,mV=.1,G=.001,Tref1; numeric degeneracy 4×4 is **not** universal SM helicity/flavor count. Grid coarse: cut12, radial12×4, nμ6, log CM12, nφ12; mid: cut12, radial12×5, nμ8, log CM16, nφ16; fine selected: cut16, radial16×5,nμ12,log CM32,nφ24. These are simultaneously changed axes, not an isolated convergence certificate.

[Executable mass/shape runner](../../../../kinetics/chi_lepton_joint_mass_gate.py); [three tests](../../../../tests/test_chi_lepton_joint_mass_gate.py); [full report](../../../../kinetics/CHI-LEPTON-JOINT-MASS-5D-AUDIT.md); [machine data](../../../../kinetics/results/chi-lepton-joint-mass-5d-2026-11-22.json).

Run: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_chi_lepton_joint_mass_gate.py; python -m kinetics.chi_lepton_joint_mass_gate --quick. Three focused tests passed locally; no RNG or copied external source. Full eight-case fine grid is an optional expensive diagnostic and is NOT required in standard CI.

## Actual positive and negative results

At mid, smooth χ1 M=.1 massless Q=8.995854126e−14, mℓ=.5 Q=8.144078344e−14 (−9.47%). Smooth χ2 M=.02 massless Q=9.066972505e−14, mℓ=.5 Q=8.185587226e−14 (−9.72%). For oscillatory synthetic f, χ1 M=.1,mℓ=.5 Qcoarse=1.64349e−15, Qmid=1.74428e−15; χ2 M=.02,mℓ=.5 Qcoarse=1.88347e−15, Qmid=1.93652e−15 and Qfine=1.97026e−15. χ2 massive-lepton oscillatory mid→fine +1.74%, **not** subpercent converged; changes multiple numeric axes and cutoff so source of discrepancy not isolated. Selected χ1 massless smooth/osc fine exists in JSON. All tested event H positive, pairwise Qχ+Qℓ=0, event energy ~floating. Common-T FD μχ=.2,μℓ=−.1 remains stationary for either χ mass and massive lepton.

## Scoped interpretation and next stage

This is a **weak functional** with two prescribed χ masses, not a time-coupled χ1+χ2+SM model. No physical χ2 abundance, φ Bose source, Zprime population, charged/neutral SM lepton and antiparticle multiplicities, SM ν/e±γ heat/entropy or FRW/H. Even conditional lepton mass effects cannot be promoted to real cosmological electron drag. No physical Γel, 86eV ΩDM, BEC, frozen f(q) or CLASS-ready. Next P0: construct TWO χ species with ONE shared evolving finite lepton bath and Pauli/energy/entropy invariants using a physically anchored weak/strong event measure, test old numerical alias, then add real φ Bose and true finite SM/gauge/FRW.

## Confirmed CI (append-only)

- [Scientific verification #45038184612](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45038184612) passed for previous Mχ1/Mχ2 finite-mℓ conditional weak source commit b08c4cae4da2d8d6cac362c13719cceee8791bb4. This does not validate physical χ2 production.
