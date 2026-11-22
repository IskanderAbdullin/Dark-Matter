# 2026-11-22 — φ Bose physical moment and inverse-fugacity Galerkin pilot

**Stage:** DM-KINETICS-20261122-PHI-BOSE-MOMENT-GALERKIN. **Run UTC beginning:** 2026-11-22T21:20:21Z measured by the isolated container clock (not inherited from earlier run). **Start branch/HEAD:** research/kinetic-equilibrium-audit-20261121 at [a6c2e975](https://github.com/IskanderAbdullin/Dark-Matter/commit/a6c2e97593bc2ea0a909b5e403b5dd406c23109c), draft PR #1; old false-thermal-boson Scientific [#45046227473](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45046227473) success. Final stage commit in Git history. Scientific program [GOAL](../../../../GOAL.md), [PLAN](../../../../PLAN.md), [kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md), [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md), [checkpoint](../../../../kinetics/research-control/STATE.json).

## Hypothesis, method, independent evidence

Previous four-pop Pauli/Bose F/R and H gates **do not** guarantee correct Bose IR particle number; 86-eV false rest-energy first-cell Nφ exceeds the real excited thermal number by ×87199/1363/22 at N16/64/256. Do not insert first-cell average as rest-energy nodal occupation: this may fix thermal N while violating event equilibrium. Proposed independent entropy-dual physical-momentum Galerkin basis with hats h_i(E), inverse fugacity η(E)=Σh_iη_i and moments M_i=∫g p²dp h_i n(E)/(2π²). Mathematically exact Σh_i=1, ΣE_i h_i=E, dS=Ση_i dM_i and symmetric negative ∂M_i/∂η_j. Use **same** h for on-shell Yukawa χ1 ↔ χ2+φ Pauli/Bose source deposition: exact dark-charge/total-energy and Jlog(F/R)≥0. This is a physical-measure instantaneous **pilot**, not a valid finite-step time integrator or a physical cosmological rate.

[Full analytical audit](../../../../kinetics/PHI-BOSE-MOMENT-GALERKIN-AUDIT.md) · [source](../../../../kinetics/phi_bose_moment_galerkin.py) · [eight focused tests](../../../../tests/test_phi_bose_moment_galerkin.py) · [compact JSON](../../../../kinetics/results/phi-bose-moment-galerkin-2026-11-22.json).

Eight scoped local tests passed in ≈0.4s; commands: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_moment_galerkin.py; OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_moment_galerkin. Deterministic NumPy/SciPy (no RNG). New GitHub CI after push requires separate confirmation.

**Positive:** at synthetic T=1GeV, mφ=86e−9 GeV, Emax12.1GeV, N16/64/256 give Nφ=0.12174521280791657/0.12174521280791631/0.12174521280791605, independently integrated Nφ=0.12174521280791706 (≲9×10⁻¹⁵ relative). Direct massless p integral 0.12179382823357307 matches ζ(3)/π². On-shell Yukawa μ1=μ2+μφ fixed point Rφ~0, chemical perturbed net φ destruction −1.162509e−5 with positive H+1.912045e−6, shell energy/charges at roundoff. Normalization of those artificial positive event weights is NOT physical.

**Required negative:** replacing original φ rest-node occupation with the real first-cell average bbar≈3.49054 instead of true b(E=m)≈1.16×10⁷ **breaks the common-T Yukawa fixed point**, giving artificial Rφ≈+0.739306 under four chosen unit-weight events, while consistent Galerkin Rφ≈0. This is not a cosmological rate. Previous false-BEC results remain permanent negative evidence.

## What remains and next executable action

Moment inversion M↦η, positivity and Bose critical boundary under **finite dt**, physically normalized Yukawa width, separately genuine condensate δ³(p), shared χ1/χ2 SM lepton/Z′ finite plasma, FRW and independent q/cutoff/time convergence are unfinished. Do not claim physical Γelastic, Ω/BEC or CLASS-ready.

**NEXT_ACTION:** implement positivity/backtracking entropy-compatible implicit Euler for a small exact on-shell Yukawa event graph by solving M(ηnew)−M(ηold)−dt R(ηnew)=0, verify dN and dU and finite-step entropy (concavity), stiff/IR and three+ dt checks, without replacing the current coupled numerical model until independently validated.
