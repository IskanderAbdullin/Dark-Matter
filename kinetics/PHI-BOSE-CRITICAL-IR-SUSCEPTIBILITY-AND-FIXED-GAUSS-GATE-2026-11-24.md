# Massive 86-eV real Bose nearcritical physical susceptibility: an underresolved Jacobian and a positive fixed-quadrature remedy (2026-11-24 UTC)

## Physical and numerical question

For χ2=.02GeV, mφ=86e−9GeV and hard Ppair=2GeV, the FIRST scalar energy interval from mφ to the exactly derived χ2+antiχ2→2φ on-shell support lower floor K=9.999751861499219e−5GeV (99.9975186keV) is collision-INACTIVE in the restricted vacuum four-channel Yukawa event network, but can contain a physically populated massive Bose gas. Future SM/FRW or additional collisions may change this support. Its thermodynamic *susceptibility* still matters independently of the restricted collision source.

Choose subcritical μφ=mφ−δ, δ>0 and T=1GeV, and the first quadratic Bernstein hat B0(E)=(1−(E−mφ)/(K−mφ))². In actual physical massive DOS, M0=∫p²dp/(2π²) B0 b and its inverse-η susceptibility G00=∫p²dp/(2π²) B0² b(1+b). Use stable E−μφ=δ+p²/(sqrt(p²+mφ²)+mφ), not subtractive E−mφ.

For δ≪mφ and p≈pc=sqrt(2mφδ), b≈T/(δ+p²/2mφ). Hence the true massive excited-state susceptibility has the positive asymptotic singularity

  G00 ~ T²(2mφ)^(3/2) / [8π sqrt(δ)]  (δ→0+, gφ=1).

This follows by p=pc x and ∫_0^∞ x²/(1+x²)²dx=π/4. The excited-state NUMBER M0 stays finite while susceptibility diverges. This is an ideal-Bose continuum mathematical critical response, NOT an actual condensate occupation or proof of a physical early-universe Bose-star population.

## Independent numerical negative and convergence

Use the same physical p²dp and Bernstein B0 as the full dark kinetic code. Compare ordinary p-Gauss n=18/36/72/144/288/576 with independent adaptive SciPy integration under u=asinh[p/sqrt(2mφδ)], and independent 50-digit mpmath in u for δ=1e−16GeV.

| δ (GeV) | Reference G00 (GeV³) | ordinary n18 G00 / reference | ordinary n18 M0 / reference |
|---:|---:|---:|---:|
| 8.6e−8 (μ=0) | 1.009868178e−6 | 1.0043379 | 1.00000151 |
| 1e−10 | 1.339831674e−6 | 0.794793 | 0.99999563 |
| 1e−12 | 3.895412784e−6 | 0.2733885 | 0.99999445 |
| 1e−14 | 2.943990148e−5 | 0.0361741 | 0.99999433 |
| **1e−16** | **2.848836855e−4** | **0.00373823** | **0.99999431** |
| **1e−18** | **2.839321416e−3** | **0.000375076** | **0.99999431** |

At δ=1e−16GeV, the physically normalized 18-point Gauss rule recovers only **0.374% of the true G00** (understates it by ≈267.5×), although the number moment seems converged to 6×10^−6. Even 576 ordinary Gauss nodes yield only **3.87% of G00**. The independent high-precision reference is

  G00 = 0.0002848836855028987709311460285161055630817 GeV³ (50-digit mpmath).

For δ=1e−18GeV, 576 ordinary nodes capture ~0.388% of G00; decreasing δ by factor100 raises the exact susceptibility by ≈9.967, agreeing with the δ^−1/2 prediction. Exact physical phase-space **capacity** and number moment convergence are not sufficient to certify nearcritical quantum Jacobian convergence.

## Constructive remedy retaining discrete entropy/moment consistency

A static **positive composite momentum Gauss** partitions p∈[0,pmax] into [0,p_min] and **20 geometrically spaced momentum intervals** between p_min=1e−16GeV and pmax=sqrt(K²−mφ²); apply n=8 Gauss on each, yielding **168 nodes**, each with positive integration weight. This node/weight grid is determined by the physical cell geometry plus declared p_min, NOT by μφ or the evolving occupancies. Using the SAME weights for the physical number/energy moments, quantum entropy and susceptibility retains the exact represented discrete derivative and KL consistency (whereas silently changing integration weights at each η without their derivatives would not).

All six δ in [1e−18,mφ] have |G00composite/G00reference−1| <**2e−7**, |M0composite/M0reference−1|<3e−10, and integrated p²dp capacity error <3e−14 relative. The deliberately coarser 12 geometric intervals ×6 Gauss nodes preserves exact p² capacity but misestimates G00 near δ=1e−16 by ~0.059%, an additional negative convergence check. The fixed p_min/resolution is only validated for these δ; smaller δ and physical condensate delta functions require new gates, not extrapolation.

The source and regression tests are:
- kinetics/phi_bose_critical_ir_susceptibility.py (stable physical-p, u-reference and fixed geometric quadrature).
- tests/test_phi_bose_critical_ir_susceptibility.py (seven tests, including 50-digit mpmath reference).
- Focused CI run #45423112252 at commit 9dfa1ff1 **success**; whole Scientific CI status independently check.

**Important:** the IR quadrature remedy is so far an independently verified *available numerical building block*, NOT yet plugged into all five coupled species. The hard-P2 first cell is legitimately collision-frozen in the limited four-channel vacuum model. Its divergent susceptibility does not by itself identify a Bose condensate, physical singular zero-mode, or change in actual cosmological relic Ω.

NEXT_ACTION: apply the geometric positive rule selectively to physical Bose IR capacity/moment/entropy within the five-field solver, compare independent q/η/thermal robustness and ensure the same quadrature is used for KL and Jacobian; test fully coupled nonthermal strong three-grid cutoffs/time, then derive physical χ+SM thermal Γ(T), FRW and actual frozen fφ(q) before CLASS likelihood.
