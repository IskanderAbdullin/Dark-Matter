# Stage 2026-11-22 — first physical Bose cell can simulate a false 86-eV condensate

- **ID:** DM-KINETICS-20261122-PHI-FIRST-ENERGY-CELL-FAKE-BEC
- **Active run UTC:** started 2026-11-22T16:22:01Z, this is a subsequent stage in the same continuous active invocation; stage-specific elapsed time not fabricated.
- **Parent research SHA:** [97224d1181c999111cfc9255ad64315bc0353681](https://github.com/IskanderAbdullin/Dark-Matter/commit/97224d1181c999111cfc9255ad64315bc0353681), draft PR #1, public Dark-Matter. Its Scientific [#45045253260](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45045253260) was in progress at last check. Previous cutoff/mesh [#45044674246](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45044674246) passed.
- **Context:** [GOAL](../../../../GOAL.md) → [AGENTS](../../../../AGENTS.md) → [kinetics specification](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md) → [invariants](../../../../kinetics/RESEARCH_GUARDRAILS.md) → [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Predeclared scientific question

Does the current real φ initial BE scalar sampled at *the rest-energy node E=mφ*, multiplied by an exact finite first energy-cell capacity, reflect the physical boson excited-state population when the mass becomes extremely light, especially the 86eV benchmark? Derive exact first-cell and full massless Bose integrals, independently compare against radial p²dp integrals for masses .1,.01,.001,.0001,86e−9 GeV and N16/64/256 at synthetic T=1GeV/Emax12.1GeV. Do not infer BEC from any b(qmin), do not change the actual quantum Yukawa solver without a separate new H/energy derivation.

## Formula and primary external context

The real scalar radial number is n=gφ/(2π²)∫_0^pmax dp p²/[exp(sqrt(p²+m²)/T)−1], with gφ=1 and μφ=0. First node b(m) at whole physical first-cell capacity w0≈g h³/(48π²) gives erroneous n0,node≈gTh³/(48π²m), while true integrated n0,true≈gTh²/(16π²): ratio ≈h/(3m)≫1 if m≪h≪T. At massless m=0, true total =gζ(3)T³/π² finite; node occupancy b(0) infinite is NOT a delta-momentum condensate. Distinct explicit condensate/normal phase-space decomposition is studied in [Blaizot et al. 2015](https://arxiv.org/abs/1510.04552). This project's count uses exactly its original density_cells energy measure; only the thermal first-cell integral is calculated independently in physical p.

## Reproducible tests and measured negative result

[Full analytic + numeric report](../../../../kinetics/PHI-IR-FIRST-CELL-ARTIFACT-AUDIT.md); [source](../../../../kinetics/phi_ir_first_cell_artifact.py); [five tests](../../../../tests/test_phi_ir_first_cell_artifact.py); [15-case compact machine data](../../../../kinetics/results/phi-ir-first-cell-artifact-2026-11-22.json).

Commands:
    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_ir_first_cell_artifact.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_ir_first_cell_artifact --full

Five local tests passed in 0.60 s, deterministic scipy.integrate.quad epsabs2e−13/epsrel2e−11. g1,T1, Emax12.1 in GeV, *not physical FRW history*. Massless full p-integral 0.12179382823357311 versus analytic ζ3/π² 0.12179382823357308. For illustrative mφ=86eV=8.6e−8GeV and N64 first true thermal cell n0,true=0.000219307612958GeV³ but naive first node n0,node=165.873845197917GeV³, a factor **756,352** too large. The full naive gas number is **1363.47×** its correct equilibrium integral. At N16/256 the total factors are 87,199×/22.289×, respectively. For the earlier toy mφ=.001GeV with N64 first-cell ratio65.75, full gas factor1.12237.

Independent first-cell-number *diagnostic only* (subtract w0 b(m), add exact p-integral of that cell) yields 86eV N16/64/256 **total ratios** 1.027080/1.003399/1.000320, a dramatic improvement. It must NOT be substituted as a b(E=m) nodal occupation into nonlinear Yukawa inverse-fugacity interpolation without preserving detailed balance and the quantum H theorem.

## Negative scope, scientific decision, next step

The finite volume **measure is correct**, but applying singular BE occupation at a single rest-energy point to the whole positive-volume cell is a provable sampling error. This is a strict falsification of using the existing coarse **thermal φ number/IR peak** as evidence for cosmic Ωφ or BEC when mφ≪h, not a no-go theorem for physical Bose dark matter or 86eV bosons. Actual generated fφ(t,q), conserved charges in full plasma, physical SM/Zprime couplings, FRW evolution, independent qmin and zero-mode scaling remain unvalidated.

Next P0: experimentally demonstrate that existing equilibrium/detailed-balance tests can pass despite this extreme initial Bose number bias; then design a correct IR-resolved, zero-mode-aware scalar representation maintaining on-shell χ1↔χ2+φ Pauli/Bose, physical radial+energy moments and quantum entropy, with at least three mesh/mass ratios. Do not add a silent Bose floor, do not identify b(m) as condensate number. Preserve negative results and checkpoint.
