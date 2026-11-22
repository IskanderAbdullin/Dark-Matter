# Stage 2026-11-22 — χ1χ2–lepton 5D phase-space cutoff independently separated from nonlinear Fermi energy-grid alias

- **ID:** DM-KINETICS-20261122-FOUR-SPECIES-CUTOFF-MESH-ORTHOGONAL
- **Run clock:** actual UTC start 2026-11-22T16:22:01Z, checked via isolated container clock. This is one stage of that active launch; exact stage duration not fabricated.
- **Branch / prior HEAD:** [74819fa4398582c2e40443772734144033bb5b0d](https://github.com/IskanderAbdullin/Dark-Matter/commit/74819fa4398582c2e40443772734144033bb5b0d), draft PR #1. Previous real [Scientific CI #45041474724](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45041474724) succeeded. New stage CI pending.
- **Project route:** [GOAL](../../../../GOAL.md) → [AGENTS](../../../../AGENTS.md) → [kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md) → [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md) → [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Question / alternatives / methods

The previous cut6→10 χ1 elastic heat difference is a combined phase-space **and mesh** change: the common energy-grid endpoint changes with elastic_cut. Distinct alternative hypotheses: (a) true finite incoming momentum-tail effect; (b) 5D weak quadrature error; (c) nonlinear energy-grid Fermi-logit interpolation alias. Previously one-axis cutoff tests could not isolate these. Introduce independent physical elastic_integration_cut ≤ elastic_cut; keep shared fermion/Bose grid, density-of-states weights, mode distributions and original FD+Yukawa source identical when changing only incoming integration region. Verify analytical FD affine-logit invariance and a physically legal nonthermal spectral counterexample.

## Reproducibility and quantitative results

Full [scientific derivation/results](../../../../kinetics/FOUR-SPECIES-CUTOFF-MESH-SEPARATION-AUDIT.md), [full source](../../../../kinetics/four_species_cutoff_mesh_separation.py), [prior 5D single-axis runner](../../../../kinetics/four_species_rate_error_budget.py), [four new targeted regression tests](../../../../tests/test_four_species_cutoff_mesh_separation.py), [two companion controls](../../../../tests/test_four_species_rate_error_budget.py) and [compact JSON](../../../../kinetics/results/four-species-cutoff-mesh-separation-2026-11-22.json).

Commands:
    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_four_species_cutoff_mesh_separation.py
    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_four_species_rate_error_budget.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.four_species_cutoff_mesh_separation --full

**Local tests**: four new tests passed in 26.11s, two companion tests passed in 13.30s; additionally nine preexisting coupled source/regressions passed after adding optional cutoff without changing defaults. No RNG, no user PC, no private CI. Prior changes in project under public PR only. GitHub CI of new commit not yet observed.

Physical synthetic case Mχ1=.1, Mχ2=.02, mφ=.001, lepton m=0, vector G=.01/mV=.1, cutoff12, incoming segmented12×Gauss3, nμ6, CM log-angle12, azimuth12; one common leptonic bath, Yukawa χ1↔χ2+real φ unchanged. Synthetic nonthermal fermion logits include amplitudes .3sin4E, −.2sin3E and +.25cos4E. Strict quantum Pauli/Bose event F/R with a real φ, not a physical relic yield.

At fixed incoming cutoff12/bins64, a **moving mesh endpoint** yields χ1 Q=−1.144123652e−11, Sdot=+8.311557324e−11, whereas a **fixed mesh cutoff20** yields χ1 Q=−1.241946092e−11 and Sdot=+5.707143000e−11. Thus the earlier mixed scan attributes ~7.9% χ heat and ~45.6% H difference to the mesh change alone, NOT physical UV tail. At cutoff8 moving vs fixed H differs ~+71.9%. A thermal FD with identical incoming phase space is *mesh-invariant* to <2e−11 relative as expected from affine logits — a strong analytic positive control, not a blanket grid artifact assertion.

With **both mesh endpoint and incoming quadrature fixed** at cut20/cut12 respectively, project 4 separate fermion energy grids N32/64/128/256 against an independently coded *analytic occupation* on-shell 5D weak integral, Qχ1=−1.267229228e−11, Qχ2=+8.750602448e−12, Hsum=+1.123620311e−10. Errors in *true* collision entropy are −86.25%, −49.21%, −18.41%, −5.73%; Qχ1 at N32 is wrong by +76.57%, N256 by −0.503%. Every node/event H identity and N/E remains at ≲5e−23 absolute. Hence exact discrete conservation/entropy does not certify actual nonthermal source rate.

## Boundary, negative archive and next gate

This is a **numerical alias and experiment-design counterexample**, not an exclusion of Bose dark matter or a physical Γ_el(T). Continuous-function "reference" uses the *same low-order 5D incoming/angle quadrature*, independent occupancy evaluation only; its continuum quadrature itself is not fully certified. Real lepton spin/gauge charge/ν channels, thermal masses, actual χ2/fφ source, finite SM e±γν/Z′ and FRW are missing. The 86-eV benchmark, ΩDM, BEC, frozen f(q), CLASS-ready remain unvalidated.

Next: independently refine Fermi mode bins and 5D incoming nlocal, nμ, CM, φ angular, physical UV cut and lepton mass; build a *separate* multi-axis error envelope for true nonthermal moments and scalar IR/thermal-source boundary. Restore full combined four-species evolving f only after numerical error gates. Maintain executable NEXT_ACTION/STATE/EVIDENCE/BACKLOG and normal commit/push, no force/merge.
