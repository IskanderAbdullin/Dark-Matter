# Stage 2026-11-22 — four coupled χ1,χ2,φ and ONE shared lepton populations

- **ID:** DM-KINETICS-20261122-FOUR-SPECIES-YUKAWA-SHARED-LEPTON
- **Active START_TIME_UTC:** 2026-11-22T15:38:09Z measured via real UTC container clock. Per-stage exact elapsed duration not separately recorded; never infer from number of tests.
- **Previous HEAD:** [4a5bb336307503631d1bde9b39d90d21353aa46e](https://github.com/IskanderAbdullin/Dark-Matter/commit/4a5bb336307503631d1bde9b39d90d21353aa46e), branch research/kinetic-equilibrium-audit-20261121 and draft PR #1; Scientific verification of this parent may still be in progress, check actual run before declaring success.
- **Context:** [root GOAL](../../../../GOAL.md), [PLAN](../../../../PLAN.md), [AGENTS](../../../../AGENTS.md), [kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md), [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md), [STATE](../../../../kinetics/research-control/STATE.json) and [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Hypothesis and competing explanations

Earlier **two χ species + one lepton** finite-vector Pauli strong operator evolved three fermionic populations, with χ2 abundance supplied synthetically. Separately the physical moving χ1↔χ2+real φ Yukawa Pauli/Bose source now evolves χ1,χ2,φ with physical radial densities. Coupling them is nontrivial: if the two generators use distinct χ grids/physical cell capacities or duplicate the lepton distribution, energy and quantum entropy are not guaranteed. Test that the **SAME four current populations** enter both collision families without hidden thermostat, that all four states evolve under the sum, and that the *correct* combined nullspace is N1+N2, N1+Nφ, Nlepton and total E, not N1 or N2 separately.

## Mathematical method and reproducibility

Both collision families are on shell: finite-vector t-channel χ1ℓ→χ1ℓ and χ2ℓ→χ2ℓ with lepton Pauli blocking, CM log-transfer and one shared lepton fℓ; physical real scalar Yukawa χ1↔χ2+φ with moving-parent Lorentz daughters, daughter Pauli, φ Bose stimulation and inverse decay. The χ1 and χ2 energy grids and w=g/(6π²)∫p²dp cell capacities are forced to match **bitwise** between the two source modules. φ and ℓ each have their own physical rest-energy lower boundary and cell measure; lepton elastic truncated incoming cutoff is 6 so the shared energy maximum12.1 safely contains Yukawa parent pmax12 and all daughter modes. Each family's event F/R reconstruction uses matching fermion logits and Bose η log inverse fugacity, paired deposition enforces exact first-energy and number moment, semidiscrete Sdot_global=Σevent W(F−R)ln(F/R)≥0. Nonlinear finite step = adaptive 2 half-Euler with full four-species error and strict quantum bounds and actual entropy; method remains FIRST order.

[New source](../../../../kinetics/four_species_yukawa_shared_lepton.py), [six focused tests](../../../../tests/test_four_species_yukawa_shared_lepton.py), [scientific audit](../../../../kinetics/FOUR-SPECIES-YUKAWA-SHARED-LEPTON-AUDIT.md), [numerical data](../../../../kinetics/results/four-species-yukawa-shared-lepton-2026-11-22.json).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_four_species_yukawa_shared_lepton.py tests/test_yukawa_bose_pauli_bounded_time.py tests/test_chi1_chi2_shared_lepton_pauli.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.four_species_yukawa_shared_lepton

**Local:** 18 tests passed together in 5.11 seconds (including parent component tests); new actual GitHub CI separately verify. The synthetic 4-species benchmark is mχ1=.1,mχ2=.02,mφ=.001,mℓ=0; y=.001, vector G=.01,mV=.1 (NOT full UV/mass/SM physical calibration), χ1/χ2/ℓ numeric g=4, φ g=1; χℓ incoming nrad4,nmu4,CM ncos8,nφ8, Yukawa parent nrad12,nangle12. End time2e7 arbitrary collision units; timestep atol2e−5.

## Physical checks and measured negative uncertainty

At bins16,Tχ1=1.1 and Tχ2=Tφ=Tℓ=1, microscopic Hdot=+7.848048e−12 with separate elastic +9.13963e−13 and Yukawa +6.93409e−12. The independently calculated quantum nodal S gradient differs from event Hdot by 1.11e−25; N/E instantaneous residuals ≤8.1e−26. Elastic transfers heat **into** the single lepton at +1.005359e−11 while Yukawa creates φ at R=+2.737842e−11. Full time has 7 accepted, 2 rejected, ΔNφ=ΔNχ2≈+0.000532534, ΔEℓ≈+0.000196965 and entropy increases +0.000153475 with N1+N2,N1+Nφ,Nℓ,energy conserved within ≤1.2e−16. Separate common T=.8 µχ1=.03,µχ2=.05,µφ=−.02,µℓ=.07 stationarity passed for mℓ=0 and .5. With G=0, φ still evolves but lepton heat vanishes: a genuine interaction-coupling negative control.

Six 3-bin × hot/cold 4-population dynamic runs show BOTH net φ production and absorption with lepton heating or cooling:

| bins | Tχ1 | accepted/rejected | ΔNφ | ΔEℓ | ΔS |
|---:|---:|---:|---:|---:|---:|
|12|.9|6/2|−.000517103|−.000139352|+.000160365|
|12|1.1|5/2|+.000535060|+.000199485|+.000154115|
|16|.9|7/2|−.000513798|−.000137761|+.000159577|
|16|1.1|7/2|+.000532534|+.000196965|+.000153475|
|24|.9|11/3|−.000511960|−.000136931|+.000158894|
|24|1.1|11/3|+.000532206|+.000195872|+.000153194|

Local time tolerance hot bins16 1e−4→2e−5→5e−6 changes ΔNφ .000533342→.000532534→.000532053 and accepted4→7→15; there is no certified global error bound. **Critical inherited negative**: χℓ vector weak rate with current nrad4 can be off by >50× for oscillatory nonthermal spectra in prior independent 5D benchmarks, despite exact Pauli and H. Present smooth thermal example is only an **internal four-population structural proof of concept**, not a physical cosmological rate or 86eV relic. φ IR bin resolution/zero mode and b=0 vacuum boundary remain absent, as do SM e±γν+Zprime thermal masses/flavor, FRW H/redshift/entropy, physical ΩDM/BEC/CLASS.

## Conclusion and actionable next task

New mathematical scope: actual four-component finite-volume χ1/χ2/φ/ONE lepton nonlinear kinetics exists and passes checked microscopic entropy, conservation, detailed-balance and both interaction-family presence. **The physical dark matter problem is NOT solved**; no measured/derived 86-eV abundance. Next P0: independent physical χℓ 5D radial/log-CM angular/cutoff and Yukawa parent/daughter IR/time bin error budget **inside** the coupled 4-component solver; then gauge/SM finite reservoir+FRW and consistent cosmic masses/couplings. Save stage, update EVIDENCE/BACKLOG/NEXT_ACTION and proceed in this session if actual time budget remains. No force or main merge.

## Parent CI check

- [Scientific verification #45040188999](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45040188999) passed on prior full three-species Yukawa stage 4a5bb336307503631d1bde9b39d90d21353aa46e. The four-species new commit's CI remains separately to verify.
