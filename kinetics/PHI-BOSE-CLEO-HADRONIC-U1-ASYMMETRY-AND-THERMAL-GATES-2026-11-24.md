# CLEO III measured hadronic A' width, exact dark U1 heavy-partner asymmetry floor, and thermal-branch inconsistency (2026-11-24 UTC)

## Motivation and three independently scoped physical gates

The full real-φ Yukawa χ1,antiχ1,χ2,antiχ2,φ plus kinetically mixed common-charge vector portal has two mandatory complications not present in a φ-only freezeout: a charged stable 20-MeV χ2 partner protected by exact net Dirac U1, and a 10-GeV visible vector mediator with real muon/tau/hadronic decays. Earlier thermal fixed-g* symmetric χ2+antiχ2 toy predicted Ωχ2pair h²≈879.959 for y=.001, far above .120. That toy deliberately required a shared SM+χ+φ temperature through freezeout; this follow-up tests whether the **explicit vector portal** even supports this crucial assumption and replaces the previously unspecified hadronic ratio with an actual primary measurement, while preserving the distinction between necessary exact bounds and conditional model outputs.

## 1. Actual published hadronic R(s) at 9.996 GeV (not exact 10 GeV)

**Primary CLEO III paper**: D. Besson et al., Phys. Rev. D **76**, 072008 (2007), DOI https://doi.org/10.1103/PhysRevD.76.072008, arXiv https://arxiv.org/abs/0706.2813, Table VII (verified in the physical table image). Its inclusive **radiatively corrected continuum Born** ratio is

    R_had(s) = σ0(e+e−→hadrons)/(4πα²/(3s)).

There are seven experimental energy points, all copied verbatim into the new numerical source. For the closest point to the illustrative **mV=10.000 GeV** vector resonance:

    sqrt(s)_CLEO = 9.996 GeV,
    R_had = 3.497 ±0.004_stat ±0.064_syst_common ±0.043_syst_uncorrelated,
    σ_R,one_point = sqrt(.004²+.064²+.043²) = 0.07720751259,
    σ_R/R = 2.208%.

**The source is 4 MeV away from mV=10.000 GeV, not at the same energy.** The experimental systematic components have different correlations across energies, and the one-point quadratic error only describes the measured point. CLEO explicitly subtracted initial-state-radiation and Υ resonance contributions; further unknown local energy interpolation, Z/A' mixing, QCD final-state threshold effects and *thermal* spectral-function corrections are NOT part of this 2.21% error.

The previous illustrative assumption R_had=3.5 has now been REPLACED in the *data-driven numerical reference* by actual measured 3.497, while retaining an API that refuses interpolation if the energy mismatch exceeds a declared 0.01GeV maximum. The default data-free branch still has R_had=None and does not fabricate physical R. Physical hadronic A' partial width is Γhad = (α ε² mV/3) R_had at leading photon-like coupling (Born convention). Coupled with e,μ,τ and both χ dark widths:

    Γ(A' dark+all charged visible at 10GeV) ≈ 0.005305164626087453 GeV,
    [Γ onepoint ±σ_R] ≈ [0.005305164625899649,0.005305164626275256] GeV,

which changes total ΓV negligibly because dark decays dominate at gD=.1,ε=1e−5. The *visible* branching ratios change substantially; the independently normalized spin-3 real-vector narrow-width (NWA) χ2pair→A'→all charged visible rate at T=1GeV is

    γ_NWA,all visible = 2.237435933268337e−15 GeV4,
    one-point statistical+systematic band [2.2108216094,2.2640502571]e−15 GeV4,
    ratio to old electron-only NWA = 6.4907477246.

These are **forward-only**, not net chemical relaxation or actual SM plasma collision rates; inverse still cancels at thermal equilibrium. If a real intermediate A' is dynamically evolved, avoid double-counting its on-shell contribution in a simultaneous χχbar↔SM+SM Breit–Wigner 2→2 source. No experimentally tested vector parameter fit is asserted.

Reproduction: [primary CLEO Table VII source](cleo2007_hadronic_R_reference.py), [physical partial-width mapping](phi_bose_vector_portal_visible_width_gate.py), [four independent regressions](../tests/test_cleo2007_hadronic_R_vector_gate.py).

## 2. A strict nonthermal, non-FRW algebraic lower bound from conserved dark U1

Under the explicitly stated minimal real-φ Yukawa and equal vectorlike χ1/χ2 U1 charges, the **seven** allowed dark Yukawa decay, Dirac/anti pair annihilation into φφ and SM, and elastic dark/e channels all preserve the integer Noether charge

    Q_D=Nχ1−Nantiχ1+Nχ2−Nantiχ2 .

φ, A' and SM carry zero dark charge. At late times the absolute net charge |QD| cannot disappear even if **all charge-symmetric pairs annihilate perfectly**. By triangle inequality the surviving number of χ1/antiχ1/χ2/antiχ2 is at least |QD|, and since the lightest charged mass is mχ2=.020GeV,

    ρ_dark,charged(today) >= mχ2 * |Q_D,comoving| / a_today³.

Normalize to TODAY'S actual SM entropy, with ηD(today)=Q_D,comoving/S_SM,comoving(today), so NO assumption that the early SM entropy separately conserved is necessary. The strict abundance floor is

    Ω_charged h² >= 0.020 GeV × [s_SM,0/(ρcrit/h²)] × |ηD(today)|,
    s_SM,0=2891.2cm^-3, ρcrit/h²=1.05375e−5 GeV cm^-3.

Requiring even just Ω_charged h²≤0.120 gives the rigorous **necessary** upper bound

    |ηD(today)| <= 2.1868082457111234e−8.

If the intended 86-eV φ occupies ≥99% of observed DM and the charged stable partner is allowed at most 1%, the necessary bound tightens to

    |ηD(today)| <= 2.1868082457111234e−10.

Unlike the earlier 879.959 toy result, this inequality is **independent of a thermal MB distribution, φ chemical potential, Bose stimulation, exact freezeout T, dark plasma kinetic equilibrium or FRW entropy-production history**, assuming only the exact U1, minimal charges/masses and chosen today's observational density bound. No measured dark-sector asymmetry ηD is supplied; the inequality does not exclude ηD=0 or validate the model. Any new genuine dark U1-violating coupling changes the model and must be separately audited. See [physical charge-source](phi_bose_dark_U1_asymmetry_floor.py), [independent tests](../tests/test_phi_bose_dark_U1_asymmetry_floor.py).

## 3. The declared electron mediator is much too weak to support shared SM thermal T at χ2 late freezeout

Independently evaluate the *exact full relativistic vacuum Yukawa χ2antiχ2→φφ* MB thermal cross-section <σv>GG (with both t/u, identity factor), multiplied by the equilibrium antiχ2 physical spin2 number, against actual RD H(T). Separately evaluate the physical **finite-mV10** t-channel χ2+e± **transport-weighted scattering-frequency proxy** νT,e/H using a different Bessel integral and the same gD=.1,epsilon=1e−5. This is NOT the true Pauli/drag eigengap and the e bath is massless MB with externally fixed g*ρ=10.75; near MeV real electron mass and neutrino/γ entropic evolution need separate treatment.

With mχ2=.02GeV:

| x=m/T | T (MeV) | exact χ2pair forward Γ/H | finite-A' electron νT,e/H (NOT drag) |
|---:|---:|---:|---:|
| 3 | 6.667 | 1060.3173 | 2.92631594e−6 |
| 5 | 4.0 | 120.7638 | 6.37710188e−7 |
| 8 | 2.5 | 4.969233 | 1.55655651e−7 |
| 10 | 2.0 | 0.610205 | 7.95767656e−8 |
| 12 | 1.667 | 0.0761003 | 4.59999518e−8 |
| 15 | 1.333 | 0.00342046 | 2.35370789e−8 |

At x≈8–10 actual dark pair chemical annihilation crosses H, yet the electron vector transport-weighted collision proxy is 6–7 orders below H. Off-shell A' χ2pair→e+e− forward rate/H is also tiny (<1e−8 at x8). The toy assumption of a *single shared SM, χ2 and φ temperature* through chemical freezeout is therefore **NOT justified by the very portal specified in the toy microphysics**. Strong χ2-φ internal thermalization or prior decoupling followed by dark entropy transfer may still yield a different physically self-consistent dark temperature; must actually be evolved. The previous conditional Ω≈880 thus remains a useful thermal-branch falsification **but is not yet a prediction of this 10GeV vector model**.

A one-parameter root of this electron νT,e/H=1 proxy with g*ρ=61.75 fixed can only estimate a mass/temperature *crossing*, NOT a model-independent kinetic decoupling temperature. The ε required to make the proxy alone equal H at x8 is enormous compared with ε=1e−5; no experimental permission is inferred, and other SM channels can alter the true thermal history. See [independent paired rates](phi_bose_chi2_thermal_branch_consistency.py), [reproducibility tests](../tests/test_phi_bose_chi2_thermal_branch_consistency.py).

## Scientific next action and limitations

The highest priority is **self-consistent χ1,antiχ1,χ2,antiχ2,realφ+SM/A' FD/BE strong kinetics plus finite SM heat and FRW**, including actual separate dark temperature if decoupled, full g*S/g*ρ, all dark/SM inverse processes and possible χ2 asymmetry. Then calculate **Ω_charged+Ωφ** and frozen PSD with independent q,dt,angle/UV and nonzero physical zero-mode checks, and only then claim any CLASS/observational verdict. At 10GeV mV, BaBar 2017 invisible single-photon study https://doi.org/10.1103/PhysRevLett.119.131804 only reports direct upper limits for mA'≤8GeV, so it **does not directly exclude or permit** this 10GeV benchmark; a correct applicable collider/beam bound remains to be checked.

References and data stored in permanent Git, not solely in ephemeral CI logs. All numerical statements above apply to the explicitly declared benchmark and scope. **class_ready=false; no physical 86eV Ω, macroscopic BEC or galaxy-halo observational model is established.**
