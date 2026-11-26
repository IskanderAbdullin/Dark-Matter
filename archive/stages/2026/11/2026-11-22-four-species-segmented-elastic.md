# Stage 2026-11-22 — full four-population χ1χ2φℓ radial/cutoff 5D elastic error gate

- **ID:** DM-KINETICS-20261122-FOUR-POPULATION-ELASTIC-SEGMENTED-P
- **Session START_TIME_UTC:** 2026-11-22T15:38:09Z, actual UTC clock. Per-stage stop clock not separately measured.
- **Branch/prior HEAD:** public IskanderAbdullin/Dark-Matter, research/kinetic-equilibrium-audit-20261121, [3e8b1069cee9d92b8369fbc9cf0cedf063c05161](https://github.com/IskanderAbdullin/Dark-Matter/commit/3e8b1069cee9d92b8369fbc9cf0cedf063c05161), draft PR#1, [real Scientific CI #45040794111](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45040794111) success on earlier four-pop stage.
- **Scientific context:** [GOAL](../../../../GOAL.md), [PLAN](../../../../PLAN.md), [kinetics spec](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md), [guardrails](../../../../kinetics/RESEARCH_GUARDRAILS.md), [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Scientific question and before-result criteria

The coupled four-state χ1/χ2/realφ/ONE-lepton time source has exact global charges/energy/H but inherits **uniform four-point Gaussian** incoming 5D χ+lepton elastic integration, previously disproved as physically accurate for nonthermal f. Test whether replacing incoming p/k with physically energy-guided segmented p Gauss, WITHOUT changing the mathematical quantum source and old default, allows independent weak-vs-strong energy/H agreement across 3+ radial local resolutions and reveals remaining UV/cutoff errors. Preserve negative oscillatory Pauli occupations (full 0<f<1) and all valid conservation laws.

## Equations and reproducible code

Use exact p(E)=√(E²−m²) at E_m=m+jΔE, maximum ΔE=1 (toy units), then apply local Gauss-Legendre directly to each physical p interval, separately for incoming χ and lepton. Existing log-transfer CM angle handles finite-vector propagator without losing F−R Pauli; same barycentric logit interpolation and χ/ℓ paired energy deposition. Update [shared χ1χ2/one-lepton 5D source](../../../../kinetics/chi1_chi2_shared_lepton_pauli.py) to accept radial_scheme and segments, forward those exact parameters from [coupled four-state source](../../../../kinetics/four_species_yukawa_shared_lepton.py). The uniform4 old default remains unchanged, preserving historical regression. Yukawa physical moving scalar source/φ Bose untouched. [Three focused tests](../../../../tests/test_four_species_segmented_elastic_gate.py), [numeric JSON](../../../../kinetics/results/four-species-segmented-elastic-gate-2026-11-22.json), [report](../../../../kinetics/FOUR-SPECIES-SEGMENTED-ELASTIC-GATE-AUDIT.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_four_species_segmented_elastic_gate.py
    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_four_species_yukawa_shared_lepton.py tests/test_chi1_chi2_shared_lepton_pauli.py

Local results: three new tests passed ~21s, ten old parent tests passed ~4s, no MC RNG. Independent WEAK benchmark using [separate 5D program](../../../../kinetics/independent_chi_lepton_5d_quadrature.py) with the same CURRENT analytic Fermi reconstructed fχ,fℓ and physically matched Gauss/angular/cutoff agrees with projected four-state strong Q to ~1e−14 relative.

## Actual positive and negative numerical results

Smooth FD(Tχ1=1.1), FD(Tχ2=Tℓ=1), φ BE(1), Mχ1=.1,Mχ2=.02,mℓ=0,mφ=.001; vector G=.01,mV=.1 and Yukawa y=.001 arbitrary common units; elastic incoming cutoff6, nμ8, log-CM ncos16, nφ16. Segments6 × nlocal3/4/5/6 yield χ1 elastic Qχ1=(−9.913486,−9.913919,−9.914279,−9.914432)e−12; microscopic elastic H=(9.012260,9.012654,9.012981,9.013120)e−13, positive. The separate χ2 elastic source remains equilibrium, while moving Yukawa Rφ is unchanged (+2.7378421233e−11). Paired charge/energy and total H exact at same collocation.

**Negative:** keeping ΔE=1, nlocal4 and all angles fixed, increasing cutoff6/8/10 gives Qχ1=(−9.913919,−11.970179,−12.586445)e−12 and H=(9.01265,10.88198,11.44222)e−13. Cut6 *misses* ~21.2% of Q relative cut10 (or cut10 is ~26.96% larger than cut6), despite radial nlocal3/4/5 stability at cut6. Both percentages have distinct denominators; quote carefully. This is a real missing UV tail, **not** evidence of χ kinetic thermalisation or physical Γel. As incoming cutoff expands, common energy Emax changes, but this particular thermal FD initial logit is affine in E, so the single-time weak Q comparison does not get raw nonlinear interpolation alias. Not true for nonthermal mode evolution.

**Explicit negative on allowed oscillatory spectrum:** with fχ1=expit(−E+.3sin4E), fℓ=expit(−E+.3cos4E), 64 energy bins and same cutoff6/angles, uniform incoming nr4 gives χ1 Qχ1≈+1.95110070e−12, segmented6 local3/4/5 gives (+.874705,+.831056,+.823638)e−12. So gross coarse source overstates the better-resolved conditional weak transfer by **2.37×**, even with true eventwise Pauli H theorem and exact energy conservation. Segmented local3→5 still ~6% difference; no sub-percent full Q for oscillatory f claimed.

## Outcome and executable next action

Partial positive: the four-population collision source now supports a physically justified incoming p/k segment rule and independently reproducible weak/strong comparison. **Definitive negative within this cutoff:** stable radial subinterval quadrature is not a physical full transfer rate: cutoff6→10 significant and nonthermal shapes still need finer radial modes. Next gate: cutoff12/16/20 at fixed ΔE, independent nμ/CM/φ, real massive e±/ν flavor, gauge charges, physical φ daughter IR and FRW H. The user-requested full microscopic 86eV bosonic DM Ω/BEC/CLASS remains unvalidated; do not promote green GitHub CI to discovery. No force, no main merge.
