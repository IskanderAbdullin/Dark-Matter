# Two χ fermions exchanging energy with ONE evolving lepton reservoir: on-shell Pauli 2→2 pilot (2026-11-22)

## Scientific need and falsifiable criteria

In a real-neutral-φ Yukawa portal χ1↔χ2+φ, both vectorlike χ species should couple to the same leptonic sector. Prior independent full relativistic 5D on-shell finite-vector weak heat evaluated χ1 and χ2 **separately**, and did not implement a common lepton fℓ(t). Merely adding two separately integrated “bath” distributions would violate the intended one-reservoir feedback. This stage builds a **single shared-lepton**, distinct-χ1/χ2 **nonlinear semi-discrete strong source** with on-shell Pauli forward/reverse event rates, then evolves all three fermion spectra simultaneously without FD thermal reset. Separate first-principles 5D weak integrations are used to verify the two strong χ heat moments and their combined entropy at the *same numeric measure*.

Predeclared criteria: each particle-number charge Nχ1, Nχ2, Nℓ must conserve separately under pure χiℓ↔χiℓ elastic 2→2; combined χ1+χ2+ℓ energy conserved, full entropy S increases, each occupation stays 0<f<1, common-T **three independent chemical potentials** are stationary, and an adaptive time step must reject actual invalid/too-inaccurate proposals. Fail closed without clipping, post-step energy projection or fixed-temperature bath reset. This stage expressly does not assert cosmological rates.

## Distinct species but only one lepton energy grid

Use Mχ1=.1, Mχ2=.02 and mℓ=0 in arbitrary energy units, vector mediator mV=.1, G=.001, synthetic fixed Tref=1, incoming 5D Gaussian orders nr=nμ=4, outgoing CM log-transfer ncos=8 and nφ=8, cutoff12. These **coarse test orders are NOT rate-converged** (previous nonthermal 5D coarse Q error ~55×), so use them only to check the structure efficiently. The historical gχ=gℓ=4 are a numerical convention, not correct leptonic helicity/flavor counts, nor a realistic portal coupling assignment.

The maximum shared energy bound is Etop=max(Mχi)+mℓ+2cut, and a single lepton grid Eℓ,j∈[mℓ,Etop] is used for **both χ scattering channels**. Separate χ grids have physically distinct thresholds Mχi but the same upper bound. Each species has positive exact 0th radial phase-space capacities

    w_i(m,g) = g/(6π²)*[(U_i²−m²)^(3/2)−(L_i²−m²)^(3/2)].

On-shell χi(E)+ℓ(ω)→χi(E′)+ℓ(ω′), with |M|² spin trace, finite-vector propagator mapped into CM log transfer, full Fermi blocking:

    F = fχi(E) fℓ(ω) [1−fχi(E′)] [1−fℓ(ω′)]
    R = fχi(E′) fℓ(ω′) [1−fχi(E)] [1−fℓ(ω)].

For both χ species, event occupancies are reconstructed from CURRENT node Fermi logits by exactly the same positive barycentric coefficients used to deposit signed microscopic event source to each χ grid and to the ONE shared lepton grid. Then, **on every event** Nχi and Nℓ are unchanged, Eχi+Eℓ is conserved, and quantum entropy production Jln(F/R)≥0. Summing χ1 and χ2 channels in the shared lepton rate yields the correct total energy/source identity, not a sum of two independent lepton baths:

    dotNχ1=dotNχ2=dotNℓ=0,
    Qχ1+Qχ2+Qℓ=0,
    dotS=∑ events W (F−R) log(F/R)≥0.

This is a semi-discrete nonlinear Pauli strong source. It is still finite-energy-bin **collocation** and uses synthetic numerical amplitude/phase space. The time integrator takes two half explicit Euler steps only when both half/full trajectories remain strictly inside Pauli bounds, full Fermi entropy does not fall, and all-species max absolute step-doubling difference ≤atol; otherwise reject/halve. This is first-order global time accuracy at best; no uniform global error certificate or cosmological Γ/H.

## Independently reproduced results

At bins=10, with two **distinct prescribed synthetic nonthermal initial distributions** (χ1 f=expit[-E/.9+.15exp(-E/.9)], χ2 f=expit[-E/.8−.10exp(-E/.5)], ℓ f=expit[−E+.05exp(−E)]), the joint instantaneous strong physical moments are

| Rate / test | Measured synthetic value |
|---|---:|
| Qχ1 | +8.469023562689831×10⁻¹⁴ |
| Qχ2 | +1.0384403681255402×10⁻¹³ |
| Qℓ | −1.8853427243945206×10⁻¹³ |
| full direct microscopic entropy | +3.4835786974671425×10⁻¹⁴ |
| energy-node logit entropy | +3.483578697467172×10⁻¹⁴ |
| χ1, χ2, lepton individual number-rate residuals | of order 10⁻²⁹ to 10⁻²⁸ |

Two **independent** calls to the earlier standalone [relativistic 5D weak quadrature](independent_chi_lepton_5d_quadrature.py) with custom callable Fermi occupations and the SAME lepton distribution return each Qχi and combined entropy to 2×10⁻¹¹ relative against the newly independently implemented paired strong code. At common T=1 and separate μχ1/T=.1, μχ2/T=−.05, μℓ/T=−.08, the joint source is stationary even for a separate tested mℓ=.5. The true density-of-states 0th weights for χ1, χ2 and ℓ sum to the analytic g pmax³/(6π²) to ≤2×10⁻¹⁴.

A real nonlinear all-three-Fermi evolution at artificial total time 10¹⁰, starting with one oversized proposal dt=10¹⁰ and local max spectral tolerance 10⁻⁷, produced **11 accepted / 7 rejected** step proposals; maximum local step-doubling difference 9.7667766874×10⁻⁸; positive minimum accepted-step entropy gain 2.1138946126×10⁻⁵; total entropy gain +0.0003400481078217332. Final individual Nχ1,Nχ2,Nℓ drifts ≤5.56×10⁻¹⁷, and total energy drift −4.44×10⁻¹⁶. All Fermi nodes remain strictly within (0,1), with smallest χ1=2.35×10⁻¹², χ2=8.26×10⁻¹⁴ and lepton=3.42×10⁻¹¹. Four focused tests passed locally (~1.7 seconds) independently of a future GitHub Actions run.

**Important negative and physical scope:** χ2's initial population above was *invented as a test fixture*, **NOT produced** from the real χ1→χ2+φ decay. All rates use severely coarse nr4/CM8 numerical checks which the previous independent 5D audit explicitly disproved as reliable quantitative heat for arbitrary nonthermal shapes. χ1 and χ2 elastic number charges **are not conserved** once their interconversion with real φ is later added; they are conserved in THIS ISOLATED ELASTIC experiment only. The shared Fermi lepton is a numerical reservoir, not the physical SM electron/positron/neutrino/photonic species network; historical g4×g4 does not justify any physical electron/ν drag or actual Γ_el(T). No φ Bose or mediator population, finite SM weak e±γν chemistry, plasma mass, FRW/H or external observations. First-order local tolerance does not establish global time/freezeout precision.

## Reproduction, archive and next P0

[Executable physically paired χ1/χ2–shared ℓ code](chi1_chi2_shared_lepton_pauli.py) · [four independent weak/strong and time tests](../tests/test_chi1_chi2_shared_lepton_pauli.py) · [compact numeric record](results/chi1-chi2-shared-lepton-pauli-2026-11-22.json) · [dated archive](../archive/stages/2026/11/2026-11-22-chi1-chi2-shared-lepton.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_chi1_chi2_shared_lepton_pauli.py

This stage's new GitHub Scientific CI must be **independently checked** after branch fast-forward push, not inferred from four passing local tests. Preserve all earlier failures of coarse 5D Q, entropy interpolation and angular/radial/cutoff convergence.

Next P0 is an independently verified **on-shell χ1↔χ2+real φ Bose/Pauli gain/loss** with Lorentz-allowed threshold Mχ1>Mχ2+mφ, exact chemical detailed balance, χ1+χ2 total fermion number, φ creation/annihilation and total energy exchange, mapped into the physically evolving common lepton bath. At the same time, actual χi+ℓ rates require full 5D q/angle/UV convergence and proper spin/flavor weights. Only after full χ1/χ2/φ/Zprime+finite e±γν SM FRW self-consistency may one attempt Γel/H, 86eV ΩDM, frozen f(q), BEC or CLASS predictions. **Scientific status: in_progress.**
