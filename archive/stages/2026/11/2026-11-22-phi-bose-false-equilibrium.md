# Stage 2026-11-22 — four-population Yukawa equilibrium can hide false 86-eV φ abundance

- **ID:** DM-KINETICS-20261122-FOURPOP-PHI-EQUILIBRIUM-ALIAS
- **Session UTC start:** 2026-11-22T16:22:01Z; this stage follows the thermal first-cell diagnostic. Exact stage elapsed time not claimed.
- **Prior HEAD:** [56552448d70a668cbbd44de99a02f82f4dd11901](https://github.com/IskanderAbdullin/Dark-Matter/commit/56552448d70a668cbbd44de99a02f82f4dd11901), research/kinetic-equilibrium-audit-20261121, draft PR #1. CI status must be rechecked.
- **Context:** [GOAL](../../../../GOAL.md), [AGENTS](../../../../AGENTS.md), [kinetics specification](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md), [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md), [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Falsifiable question and alternate interpretation

A separate first-cell Bose quadrature proof showed an artificial divergence in total φ number, even though the scalar F/R gain/loss is exactly microreversible. Test whether a **complete shared χ1–χ2–φ–one-lepton solver** can exhibit correct chemical detailed balance and microscopic H while integrated scalar number remains grossly wrong. The alternative that equilibrium correctness certifies total density is directly falsified if both coexist.

## Methods, controls and data

At common T=1 and μ1=μ2=μφ=0 for χ1(.1),χ2(.02),real φ(m=.001 or86e−9 GeV),lepton (massless), set elastic G=0 to isolate moving Yukawa but retain the existing coupled state/ledger/Pauli-Bose 1↔2 RHS. Six independent deterministic states (2 masses × N16/64/256 energy cells) are compared against independently integrated exact radial physical thermal p²dp number with same scalar endpoint. No FRW, no nonthermal produced φ, no true condensate mode.

[New driver](../../../../kinetics/phi_bose_false_equilibrium_gate.py), [2 passing local tests](../../../../tests/test_phi_bose_false_equilibrium_gate.py), [JSON](../../../../kinetics/results/phi-bose-false-equilibrium-2026-11-22.json), [science report](../../../../kinetics/PHI-BOSE-FALSE-EQUILIBRIUM-AUDIT.md). Run OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_false_equilibrium_gate.py and python -m kinetics.phi_bose_false_equilibrium_gate.

## Positive and negative evidence

All six real coupled-model F/R tests return |Rφ|≲5e−26, entropy production≲4e−41, energy residual≲3e−41, stoichiometric residual≲2e−41, confirming true microscopic chemical fixed points on this *finite grid*. However thermal φ densities for m86e−9 T1GeV: at N16 Nφ=10616.026 versus true0.1217452 (**87199×**), N64 Nφ165.9958 versus true0.1217452 (**1363.47×**), N256 Nφ2.71357 versus true0.1217452 (**22.289×**). Toy mφ=.001,N64 overcounts 12.24%. These are synthetic GeV³-equivalent normalized model numbers, not observed cosmological abundances.

**Scientific negative:** exact event Pauli/Bose conservation, thermal detailed balance and vanishing source/entropy do NOT imply a correct *integrated scalar density* when the first energy-node occupation diverges as mφ→0. Cannot claim 86eV ΩDM, BEC or CLASS from this solver's bφ low node.

## Next reproducible action

Derive momentum/IR-cell representation with genuinely physical number and energy integrals, independent zero-condensate component if justified, then prove Pauli/Bose F/R and H theorem under the same positive transfer weights; examine exactly zero boson boundary and 3+ IR/tail/time grids; do not silently replace b(E0=mφ) by the average first-cell occupation. Check actual CI after push and update EVIDENCE/BACKLOG/NEXT_ACTION. The result is a numerical pathology, NOT a particle-physics no-go theorem.
