# First coupled χ1,χ2,real φ and ONE common lepton kinetic population with Pauli/Bose (2026-11-22)

## Research question and independent gates

Previously (1) a physically on-shell Lorentz-boosted χ1↔χ2+real φ **1↔2** operator evolved three physical-energy Fermi/Bose populations, and (2) a separate finite-vector χ1ℓ and χ2ℓ **2↔2** full-Pauli elastic pilot evolved two χ populations but only ONE shared lepton. Combining them naively can double-count lepton baths, mismatch particle capacities and destroy the global charge/energy/entropy nullspace. The first priority is a mathematically consistent **four-population nonlinear** coupled kinetic source, not a cosmological relic prediction.

Before the calculation, require (i) EXACT common χ1 and χ2 physical-energy grids and radial cell densities in BOTH collision generators, with a single lepton fℓ, (ii) separate contributions to H from all three 2→2/1↔2 channels, nonnegative total entropy, (iii) conserved Nχ1+Nχ2, Nχ1+Nφ, Nℓ and combined E at the semidiscrete level and after every accepted time step, (iv) nonzero finite lepton heat and simultaneous net boson production/absorption, (v) full common-T FD/BE chemical fixed point at independent µχ1,µχ2,µℓ,µφ subject to µχ1=µχ2+µφ, and (vi) at least three bin grids and time tolerances, preserving negative error estimates.

## Physical event structure and shared node equations

For each χ species i∈{1,2}, use the previously independently checked finite-vector t-channel 5D on-shell χi+ℓ↔χi+ℓ collision with exact Pauli forward and inverse terms and CM forward-log-transfer variable. Its lepton source from both channels is added **once into the same current lepton population**; neither channel may reset the lepton temperature or assume fixed FD after a step. The physical real-scalar Yukawa χ1↔χ2+φ forward/inverse is the previous boosted scalar 1↔2 Fermi/Bose on-shell operator, now with the **same current χ1 and χ2 occupations** as the elastic collisions. The full system is:

    ∂t fχ1 = Cχ1ℓ[fχ1,fℓ] + CY,χ1[fχ1,fχ2,bφ]
    ∂t fχ2 = Cχ2ℓ[fχ2,fℓ] + CY,χ2[fχ1,fχ2,bφ]
    ∂t fℓ  = Cℓχ1[fχ1,fℓ] + Cℓχ2[fχ2,fℓ]
    ∂t bφ  = CY,φ[fχ1,fχ2,bφ].

Each species physical energy node E_{s,j} is normalized to its clipped radial cell capacity w_{s,j}=g_s/(6π²)[(U²−m_s²)^(3/2)−(L²−m_s²)^(3/2)]. The χ1 and χ2 node grids and w arrays are **bitwise shared** between the two generators (not copied from incompatible simulation domains). Fermion logit reconstruction and boson log inverse fugacity use the identical nonnegative first-energy-exact barycentric weights as event deposition; thus the isolated 2→2 and 1↔2 quantum H theorems add.

At every instant the left nullspace of the FOUR-species coupled stoichiometry contains

    dot Nχ1 + dot Nχ2 = 0       (dark fermion charge)
    dot Nχ1 + dot Nφ  = 0       (closed Yukawa scalar charge combination)
    dot Nℓ            = 0       (elastic lepton number)
    dot Eχ1+dot Eχ2+dot Eℓ+dot Eφ = 0.

In contrast, **Nχ1 and Nχ2 individually are NOT invariant** with Yukawa conversion; erroneously enforcing the old elastic per-species checks would be a serious physical error. Combined instantaneous Sdot_nodes = Hdot_χ1ℓ + Hdot_χ2ℓ + Hdot_Yukawa = Σevent W(F−R)ln(F/R) ≥0. The time stepper is again two recomputed explicit half-Euler steps, accepting only physical Fermi (0,1), Bose >0 and actual nondecreasing entropy with a full-spectrum error estimator; this remains **FIRST ORDER** and is not a certified global time/IR method.

## Synthetic experiments and positive evidence

Shared toy energies: mχ1=.1,mχ2=.02,mφ=.001,mℓ=0, real Yukawa y=.001, vector parameters G=.01,mV=.1 (effective numerical, not UV-calibrated), Fermi spin factors historically gχ1=gχ2=gℓ=4, real φ gφ=1; lepton elastic incoming cutoff6 giving physical energy grid max12.1, moving Yukawa parent momentum cutoff12 (guaranteed inside same grid). Elastic p/k/μ Gauss4×4×4, CM polar8 and azimuth8 with log transfer, Yukawa mother Gauss12 and rest angle12. All energies and t are arbitrary *synthetic* units; numerical end time2×10⁷, initial dt10⁷, full-spectrum atol2×10⁻⁵.

With initial FD χ1 temperature1.1 and FD χ2/lepton temperature1, BE φ temperature1, bins16:
- Full instantaneous Hdot=+7.84804797009×10⁻¹², of which Yukawa +6.93408540454×10⁻¹² and **elastic +9.13962565558×10⁻¹³**.
- Independent nodal entropy gradient vs event H differs by only 1.11×10⁻²⁵; charge rate residuals ≤3.2×10⁻²⁶ and total energy rate residual ≈−8.1×10⁻²⁶.
- Yukawa initial net Rφ=+2.73784212332×10⁻¹¹; elastic heat flow from hot χ1 to the shared lepton is +1.00535882211×10⁻¹¹, with χ2 elastic initially at common T stationary.
- At time2×10⁷ the four populations were truly updated at 7 accepted and 2 rejected steps: ΔNφ=ΔNχ2≈+0.00053253438, lepton energy ΔEℓ=+0.00019696457, actual ΔS=+0.00015347492; conserved combined three numbers and energy to ≤1.2×10⁻¹⁶, all physical bounds held.

An independent full nonzero-µ fixed point at common T=.8 uses µχ1=.03, µχ2=.05, µφ=−.02, µℓ=.07. For both mℓ=0 and massive mℓ=.5, the maximum absolute node time derivative is below 6×10⁻²⁴ and Hdot~10⁻³⁹. Setting G=0 isolates the Yukawa network: φ still evolves, but shared lepton energy is unchanged to floating accuracy, confirming the difference between a genuinely coupled bath and a decorative extra lepton state.

**Six real time runs** under three energy grids and two distinct mother temperatures:

| bins | Tχ1 | accepted/rejected | net ΔNφ (also ΔNχ2) | ΔEℓ | ΔS |
|---:|---:|---:|---:|---:|---:|
| 12 | 0.9 | 6/2 | −0.000517103 | −0.000139352 | +0.000160365 |
| 12 | 1.1 | 5/2 | +0.000535060 | +0.000199485 | +0.000154115 |
| 16 | 0.9 | 7/2 | −0.000513798 | −0.000137761 | +0.000159577 |
| 16 | 1.1 | 7/2 | +0.000532534 | +0.000196965 | +0.000153475 |
| 24 | 0.9 | 11/3 | −0.000511960 | −0.000136931 | +0.000158894 |
| 24 | 1.1 | 11/3 | +0.000532206 | +0.000195872 | +0.000153194 |

Thus the actual coupled system permits both directions of real φ conversion and BOTH signs of heat flow to the **single** lepton bath. All six cases conserve the proper global nullspace and produce entropy, but no absolute cosmic time/rate follows.

A first-order temporal tolerance screen at bins16,Tχ1=1.1 (local tol10⁻⁴,2×10⁻⁵,5×10⁻⁶) changes accepted/rejected steps 4/1,7/2,15/5 and final ΔNφ +0.0005333420,+0.0005325344,+0.0005320529. Thus a single local atol and green energy/H are not numerical convergence certificates.

## Explicit negative limitations, literature and follow-up

The combined discrete dynamics is still **inadequate for physical 86-eV dark matter**. Its χℓ component inherits a particularly coarse nrad4/nμ4 on-shell 5D integration whose shortcomings were established in prior negative nonthermal tests (Q sign flips and >50× coarse alias for oscillatory spectra); our initial smooth FD parameters make this pilot interpretable as a *structural test only*. Even though bins12/16/24 show a few-percent stable net number, φ very low-energy Bose bins, exact b=0 boundary, IR zero mode and high-momentum tail have NOT passed proper continuum convergence. Lepton ν flavor/spin multiplicities, charged electron mass/plasma screening, physical Z′ mass/gauge symmetry/thermal masses, chemical e±γν processes, FRW/H(g*)/expansion/redshift and physical dimensionful couplings are absent. No full consistency with Lyman-α, BBN/ΔNeff, CMB, relic Ω or macroscopic condensate has been shown.

Related independent theory: [Yin 2023](https://arxiv.org/abs/2301.08735), [Sakurai & Yin 2025](https://arxiv.org/abs/2410.18968). The present method and tests were implemented independently; no external code copied.

Reproduce [four-species source](four_species_yukawa_shared_lepton.py), [6 focused tests](../tests/test_four_species_yukawa_shared_lepton.py), [machine data](results/four-species-yukawa-shared-lepton-2026-11-22.json), [archival stage](../archive/stages/2026/11/2026-11-22-four-species-shared-lepton-yukawa.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_four_species_yukawa_shared_lepton.py tests/test_yukawa_bose_pauli_bounded_time.py tests/test_chi1_chi2_shared_lepton_pauli.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.four_species_yukawa_shared_lepton

**Next scientific gate:** independent 3+ physical discretization controls in combined kinetics (χℓ p/k/CM log-angle and φ mother/daughter energy/IR/time/cutoff); exact zero-Bose endpoint and Bose zero-mode treatment; gauge-complete SM ν/e±γ bath, Z′ and scalar breaking sector coupled to FRW/H. Only after valid dimensional rates/physical 86eV masses and long-time freezeout/perturbations could ΩDM/BEC/frozen f(q)/CLASS-ready be considered. Status *in progress*.
