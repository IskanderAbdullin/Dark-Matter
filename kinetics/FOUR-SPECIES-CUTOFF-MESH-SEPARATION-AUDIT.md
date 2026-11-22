# Orthogonal 5D incoming cutoff versus χ1–χ2–φ–lepton energy-mesh resolution (2026-11-22)

## Research question and falsification gate

A previous coupled four-population numerical audit varied \`elastic_cut\` at constant number of energy nodes. In the implementation this parameter simultaneously controlled (i) the incoming relativistic χℓ 2→2 p,k phase-space cutoff and (ii) the upper endpoint/spacing of **every** physical χ1,χ2,φ,lepton energy mesh. Those are independent approximations. On a Fermi–Dirac spectrum the Fermi logits are affine in energy, so grid interpolation happens to be exact; on a genuinely nonthermal spectrum changing the grid can alias the physics by a large amount even when species numbers, total energy and quantum entropy are all mathematically conserved.

The prespecified experiment is to introduce a backwards-compatible separate \`elastic_integration_cut\` smaller than or equal to the shared mesh cutoff. Hold incoming radial Gauss segments (fixed maximum ΔE≈1), relative and outgoing angles, Yukawa source, physical coupling and **the same mesh** fixed while scanning the incoming cutoff. Repeat by moving the mesh with the cutoff (old behavior). At fixed incoming phase-space cutoff, compare node-interpolated strong Qχ1/Qχ2/H with a separately implemented **continuous analytic-occupation 5D weak reference**, at ≥4 energy-mesh resolutions.

## Mathematical and implementation controls

Two vectorlike Dirac species χ1 M=.1, χ2 M=.02, a real φ mφ=.001 and ONE common lepton with illustrative mass zero, all dimensionless synthetic, on-shell finite-vector χℓ scattering at G=.01,mV=.1 plus moving χ1↔χ2+φ Yukawa. In the existing finite-volume scheme fermion node F/R are reconstructed from **piecewise affine logit(f)** and deposited with the same barycentric energy coefficients. Thus exact eventwise Qχ1+Qχ2+Qℓ=0 for elastic channels, full total energy and Ndark,N1+Nφ,Nℓ are conserved, and Sdot=ΣW(F−R)ln(F/R)≥0. These structural identities do not imply the reconstructed logit spectrum is accurate.

Added \`elastic_integration_cut\` to the coupled configuration and \`incoming_cut\` to the finite-vector 5D strong RHS; if unspecified the historical integration cutoff = mesh cutoff behavior is preserved. If incoming cutoff exceeds physical energy-mesh coverage the implementation raises ValueError rather than extrapolating an energy mode. The selected physical mesh \`mesh_cut=20\` has upper E≈40.1, while incoming cuts 8,10,12,16,20 are varied at constant segmentation ΔE≈1 and fixed angular rules. There is no FD reset or post-integration entropy repair.

Two nonthermal counterexample spectra, all strict 0<f<1:
  fχ1(E)=expit[−E/1.1+0.3sin(4E)]
  fχ2(E)=expit[−E/.9−0.2sin(3E)]
  fℓ(E)=expit[−E+0.25cos(4E)].
The real φ starts with Bose equilibrium BE(Tφ=1); this is an initial synthetic source, not a cosmological yield or evidence of BEC. Incoming p,k use 12 physical energy-guided p segments with local Gauss3, nμ6, outgoing CM log-transfer polar12 and azimuth12. Reproducible and fully deterministic, no RNG or external observational data.

## Positive controls and decisive negative numbers

**Positive special case:** For common thermal FD initial χ1,χ2,lepton occupations and the exact same incoming cutoff8, using a 32-bin mesh with its upper endpoint moving from 16.1 to 40.1 changes Qχ1, Qχ2, elastic entropy and moving Yukawa Rφ only at floating precision (<2×10⁻¹¹ relative). This confirms the expected **affine Fermi-logit exact-interpolation limit** and prevents blaming all earlier FD cutoff trends on the grid artifact.

**Negative nonthermal result:** For exactly the same incoming cut12 and 64 nodes, the historical *moving* energy mesh (upper 24.1, h=.375) yields Qχ1=−1.14412365×10⁻¹¹ and elastic Sdot=+8.31155732×10⁻¹¹. Holding mesh cutoff20 fixed (upper 40.1,h=.625) gives Qχ1=−1.24194609×10⁻¹¹ and Sdot=+5.70714300×10⁻¹¹. The inferred heat differs by ~7.88%, the entropy by **+45.6%**, without changing incoming phase space, microscopic Yukawa, couplings, or angular rule. At cut8 the entropy-moving-mesh bias rises to ~+71.9%; at cut16 it remains ~+30.3%. Thus historical nonthermal cutoff trends can combine **two independent numerical errors** and should not be interpreted as physical tail truncation alone.

**Direct continuum-occupation 5D weak cross-check at exactly matched incoming quadrature (cut12):** independently evaluated Qχ1=−1.2672292278252052×10⁻¹¹, Qχ2=+8.75060244836362×10⁻¹², summed collision entropy +1.1236203110330876×10⁻¹⁰ (masses .1/.02 separately, shared analytical lepton f). At fixed mesh upper 40.1, four positive-logit interpolation energy grids are:

| bins per species | h | χ1 Q rel. error | χ2 Q rel. error | elastic Sdot rel. error |
|---:|---:|---:|---:|---:|
| 32 | 1.2500 | +76.57% | −15.76% | **−86.25%** |
| 64 | 0.6250 | −2.00% | −2.13% | **−49.21%** |
| 128 | 0.3125 | −1.96% | −0.59% | **−18.41%** |
| 256 | 0.15625 | −0.503% | −0.230% | **−5.73%** |

All these local comparisons have event-vs-node entropy-identity and total energy residuals at ~10⁻²⁴ absolute, with positive microscopic H. Even 256 modes therefore do **not** establish subpercent physical entropy convergence. A coarse 32-bin scheme can be energy/charge conserving while producing ~86% too little *true* collision entropy for the oscillatory f.

The direct analytic-occupation weak evaluation is an **independent code route at identical incoming/angle quadrature**. It does NOT certify the radial/CM/angular/cutoff quadrature itself. It only isolates logit-grid reconstruction error; no physical Γelastic, FRW or 86-eV relic rate is inferred. For physical long-time claims, independently vary energy-mode grids and phase-space quadrature, and calculate convergence of fχ1, fχ2, fφ and lepton distributions, not merely their integral charge constraints.

## Reproduction and permanent evidence

- [Modified physical χ1/χ2 shared-lepton 5D strong kernel](chi1_chi2_shared_lepton_pauli.py) with independent incoming cutoff.
- [Modified four-population coupling](four_species_yukawa_shared_lepton.py) retaining Yukawa and common lepton.
- [Single-axis 5D rate scanner](four_species_rate_error_budget.py), [matched direct-weak comparison](four_species_cutoff_mesh_separation.py).
- [Targeted tests](../tests/test_four_species_cutoff_mesh_separation.py), [complementary rate-axis regressions](../tests/test_four_species_rate_error_budget.py).
- [Compact machine data](results/four-species-cutoff-mesh-separation-2026-11-22.json), [dated stage archive](../archive/stages/2026/11/2026-11-22-four-species-cutoff-mesh-separation.md).

Run:
    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_four_species_cutoff_mesh_separation.py tests/test_four_species_rate_error_budget.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.four_species_cutoff_mesh_separation --full

The four new focused controls and two prior rate-axis checks passed **locally** in separate runs; real GitHub Scientific CI for the eventual new commit must be inspected separately. Earlier Scientific CI #45041474724 on the prior four-species segmented stage succeeded, but does not validate this new separation.

**Next P0:** hold physical mode coverage and bin resolution independently while scanning cutoff8/12/16/20 and nlocal/nμ/CM/φ angular axes and a finite electron mass; quantify conditional Q/S source error separately from time discretization. Then physical ν/e±γ/Z′ plasma, real χ masses/couplings, H/FRW, low-energy φ boundary, fφ(q), ΩDM and structure/CLASS. The 86-eV benchmark remains **unvalidated** and CLASS-ready=false.
