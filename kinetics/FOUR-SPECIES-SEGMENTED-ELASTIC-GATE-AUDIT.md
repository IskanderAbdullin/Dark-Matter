# Coupled four-population χ1χ2φℓ: remove incoming radial alias from elastic 5D Pauli collision, expose UV cutoff (2026-11-22)

## Why this is a necessary physical—not cosmetic—test

The four-population nonlinear χ1,χ2,real φ,ONE lepton simulator has conserved global N1+N2,N1+φ,Nℓ,energy and H in its *first* scoped run, but its finite-vector χℓ elastic collision was evaluated with the known-underresolved **uniform four-point p and k Gaussian integration**. Earlier independent full 5D Pauli tests established that for legal nonthermal spectra a 4-point incoming rule can change the sign of heat and overestimate the better-resolved weak energy rate by ~55×. Green conservation/entropy checks or good Yukawa physics do not erase an incorrect elastic rate. This stage integrates the independently tested energy-guided momentum segmentation into BOTH χ1ℓ and χ2ℓ channels feeding the SAME lepton—and explicitly separates radial, outgoing angle and finite cutoff.

**Predeclared scientific checks:** preserve old default for backward comparability; make incoming scheme and segment count explicit in four-population config; at identical current nonthermal Fermi arrays compare the projected strong elastic Qχ and entropy against a *separately implemented* 5D on-shell weak integral on EXACTLY the same full p,k,μ,log-CM,azimuth measure; vary ≥3 per-segment Gauss orders and ≥3 kinetic cutoffs without reducing energy-segment resolution; test an oscillatory physical Pauli counterexample against the legacy uniform4 source. Retain actual failures.

## Algorithm and independent reference

For each species of rest mass m and incoming physical kinetic cutoff Λ, choose fixed energy edges E_j=m+jΔE with ΔE≤1 in an arbitrary T=1 scale, convert them to p_j=√(E_j²−m²), and apply Gauss-Legendre with nlocal p nodes separately on each [p_j,p_{j+1}]. This has the exact ∫p²dp radial volume with nlocal≥2, preserving the physical E(p) threshold. Use this for both incoming χ and lepton radial integrals in the same 5D Lorentz-boosted finite-vector scattering. The existing CM logarithmic transfer integration removes the leading [mV²/(mV²−t)]² forward peak but does **not** solve the UV/tail problem. At every event the microscopic full F/R Pauli, positive event weights and same-logit paired χ/lepton deposition retain exact separate χ1/χ2 and lepton numbers, χ energy+lepton exchanged energy, and entropy per channel.

The new [four-species code](four_species_yukawa_shared_lepton.py) forwards the separately configurable elastic_radial_scheme and elastic_segments to [the existing shared-lepton source](chi1_chi2_shared_lepton_pauli.py), without changing any old uniform_p default or the already derived boosted Bose Yukawa. The χ1/χ2 and φ physical energy cell capacities remain shared with the moving Yukawa source; do not confuse this mathematical improvement with a calibrated physical gχ/gℓ coupling or actual cosmic plasma.

## Independent real numerical evidence

Model input *strictly synthetic*: mχ1=.1,mχ2=.02,mφ=.001,mℓ=0, G=.01,mV=.1,y=.001; initial FD χ1 T1=1.1, FD χ2/lepton T=1, φ BE T=1; χℓ incoming cutoff Λ=6, Yukawa mother pmax12, physical common energy grid covers Emax12.1. For the elastic rate test use nμ=8, mapped CM polar16, azimuth16 and energy-guided six segments; nonthermal f1 and fℓ are recomputed using the *current same logit energy node distributions*. Independently compare with [the original standalone 5D weak evaluator](independent_chi_lepton_5d_quadrature.py), not the same RHS called twice.

For three and four local orders:

| incoming six segments × nlocal | χ1 elastic Qχ1 (arbitrary energy-rate units) | elastic Hdot |
|---:|---:|---:|
| 6 × 3 | −9.91348551603×10⁻¹² | +9.01225956003×10⁻¹³ |
| 6 × 4 | −9.91391897416×10⁻¹² | +9.01265361287×10⁻¹³ |
| 6 × 5 | −9.91427924172×10⁻¹² | +9.01298112884×10⁻¹³ |
| 6 × 6 | −9.91443222265×10⁻¹² | +9.01312020241×10⁻¹³ |

Standalone independent weak Qχ1 agrees with each projected strong first-energy moment at ~1e−14 relative and positive entropy agrees to roundoff; χ2 elastic is stationary since its T matches the ONE lepton. The other moving Yukawa source keeps initial Rφ=+2.737842123319×10⁻¹¹ **unchanged**. However at the old uniform nr4 and nμ4, ncos8 and nφ8, χ1 Qχ1≈−1.00535882211×10⁻¹¹, illustrating why comparing different angular and radial measures is not a clean rate convergence test.

**Real unresolved UV negative:** retaining fixed ΔE=1, nlocal4 and nμ8/ncos16/nφ16 while moving the physical incoming cutoff 6→8→10 gives

| elastic cutoff | χ1 Qχ1 | elastic Hdot |
|---:|---:|---:|
| 6 | −9.91391897416×10⁻¹² | +9.01265361287×10⁻¹³ |
| 8 | −1.19701785040×10⁻¹¹ | +1.08819804581×10⁻¹² |
| 10 | **−1.25864446218×10⁻¹¹** | **+1.14422223835×10⁻¹²** |

Thus a numerically stable cutoff6 radial integral is still **~26.96% smaller in |Q|** than cutoff10. Cutoff6→8 alone is a large change, and 8→10 still ~5.15%. All these integrals satisfy entropy and charge/energy closure; these checks are *not* a physical Γelastic calibration. Increasing cutoff also changes the nominal Emax of the finite node grids, but the thermal FD logits are exactly affine, so this particular *initial-state weak* comparison does not introduce a spectral interpolation bias; it would be invalid to make the same claim for arbitrary nonthermal f.

**Explicit oscillatory failure control:** choose legal fχ1(E)=expit(−E+.3sin4E), fℓ(E)=expit(−E+.3cos4E), with χ2 and φ separately positive. At the same incoming cutoff6, angular nμ8/log-CM16/azimuth16 and 64 physical energy bins:
- uniform incoming nr4 gives Qχ1=+1.95110070325×10⁻¹²;
- segmented6×3 gives +8.74704911957×10⁻¹³;
- segmented6×4 gives +8.31056364988×10⁻¹³;
- segmented6×5 gives +8.23637721126×10⁻¹³.

Uniform four-point under this nonthermal state **overestimates Q by 2.37×**, although both schemes have positive full Pauli entropy and conserved moment rates. Segmentation at orders3→5 is not yet exact; the nonthermal φ spectrum itself has not been cosmologically generated. The 3 test functions dedicated to this physical-rate discrepancy passed locally in ~21s; previous 10 tests of both original coupled/elastic sources remained passing when using their unchanged default.

## What remains, and next executable gate

The stage demonstrates that forward segmented momentum measures are usable in the genuine combined on-shell χ1/χ2 source, and that a physically large cutoff error can survive even when radial + angular discretizations are locally accurate. It **does not certify** full 2→2 χ–SM weak/strong convergence or the physical Yukawa/φ IR production rate. In the next stage extend cutoff12/16/20 at fixed segment energy width, cross-check 3+ μ/CM/azimuth and both χ masses+massive electron channels. Then independently refine moving mother p and real φ daughter low-energy E bins and true bφ=0 boundary. Need complete ν_e,ν_μ,e±γ/Zprime thermal plasma and H/g*(T)/FRW and experimental constraints before 86eV φ freezeout, ΩDM/BEC or CLASS-ready can be claimed.

Archived artifacts: [modified shared-lepton kernel](chi1_chi2_shared_lepton_pauli.py), [combined four-state code](four_species_yukawa_shared_lepton.py), [3 focused tests](../tests/test_four_species_segmented_elastic_gate.py), [compact numeric results](results/four-species-segmented-elastic-gate-2026-11-22.json), [permanent stage](../archive/stages/2026/11/2026-11-22-four-species-segmented-elastic.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_four_species_segmented_elastic_gate.py
    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_four_species_yukawa_shared_lepton.py tests/test_chi1_chi2_shared_lepton_pauli.py

No physics discovery, no physical Γelastic or cosmological 86-eV mass inferred.
