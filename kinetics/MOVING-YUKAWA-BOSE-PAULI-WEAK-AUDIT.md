# Fully moving χ1↔χ2+real φ Yukawa Bose/Pauli weak source: threshold, inverse channel and thermal sign (2026-11-22)

## Question, independent physical basis and test criteria

The newly built two-χ vector Pauli elastic source evolves χ1, χ2 and **one shared lepton bath**, but χ2 had an *arbitrarily prescribed initial distribution* rather than one produced by the real φ Yukawa conversion. The prior project [microphysical model](MICROPHYSICS-MODEL.md) already derived the rest-frame scalar Yukawa decay width and checked pointwise Fermi/Bose thermal gain–loss. This stage extends those checks to a **moving mother and the fully on-shell, boosted daughter energies**. It remains a prescribed-spectrum WEAK number/heat source, NOT nonlinear mode-resolved three-species decay or a cosmological relic calculation.

Falsifiable gates: (i) physical two-body threshold and exact massless-daughter width, (ii) Lorentz boosted daughter energy and 4-momentum conservation at multiple mother momenta, (iii) correct *forward AND inverse* Pauli/Bose source with exact common-temperature chemical equilibrium at nonzero μ, (iv) species source stoichiometry dotNχ1=−dotNχ2=−dotNφ, (v) Qχ1+Qχ2+Qφ=0 and nonnegative true microscopic entropy, (vi) distinguish zero initial boson inverse from finite Bose stimulation and **preserve an actual negative physical direction** when χ1 is colder than daughter bath, (vii) at least three radial orders for one same physical synthetic fixture, not a generic numerical accuracy claim.

## Real scalar tree Yukawa trace, kinematics and quantum gain/loss

Use the renormalizable low-energy operator

    L_Y = −y φ (barχ2 χ1 + barχ1 χ2)

with a **one-component real φ** neutral under U(1)′, and vectorlike χ1/χ2 of the SAME U(1)′ charge (otherwise the original scalar Yukawa is forbidden; see [gauge-charge audit](REAL-PHI-YUKAWA-GAUGE-CHARGE-AUDIT.md)). The low-energy Z2 parities differ for χ1 and χ2. For fixed masses m1>m2+mφ and dimensionless y,

    λ = [m1²−(m2+mφ)²] [m1²−(m2−mφ)²]
    |M_Y|²_avg = y²[(m1+m2)²−mφ²]
    Γ_rest = |M_Y|²_avg sqrt(λ)/(16π m1³).

The physical threshold is exact: Γ_rest=0 for m1≤m2+mφ. At m2=mφ=0, Γ_rest=y²m1/(16π), independently checked.

For a mother with laboratory momentum p and energy E1=sqrt(m1²+p²), daughter rest-frame energies E2*=(m1²+m2²−mφ²)/(2m1), Eφ*=(m1²−m2²+mφ²)/(2m1), p*=sqrt(λ)/(2m1). With isotropic mother-rest decay angle c=cosθ*,

    E2= (E1/m1) E2* + (p/m1) p* c
    Eφ= (E1/m1) Eφ* − (p/m1) p* c
    E2+Eφ=E1.

The mother-spin-averaged differential probability is uniform dc/2 for this scalar channel; laboratory time dilation is Γ_rest m1/E1. For arbitrary prescribed isotropic fermion functions f1(E1), f2(E2) and real boson bφ(Eφ):

    F=f1(1−f2)(1+bφ), R=(1−f1)f2 bφ, J=F−R,
    dR/dpdc = [g1 p²/(2π²)] Γ_rest(m1/E1) J/2,
    dotN1=−dotR, dotN2=dotR, dotNφ=dotR,
    Q1=−∫E1 dR, Q2=∫E2 dR, Qφ=∫Eφ dR,
    Sdot=∫(dR/J) J ln(F/R) ≥0.

For zero initial φ, R=0 exactly and forward production still exists. The formal entropy derivative at strictly zero Bose occupation is **positively infinite**, not a finite value to be silently clipped; no full entropy-at-boundary solver is claimed. For common T, the correct detailed balance is μ1=μ2+μφ, and φ chemical potential need not be zero in a CLOSED pure 1↔2 network. In a complete model other scalar number-changing interactions may force μφ=0. Physical Bose condition is μφ<min Eφ, not merely a chemical affinity equation.

## Independent synthetic moving numerical checks

All dimensionless energies here have **no direct observed mass interpretation**: illustrative m1=.1,m2=.02,mφ=.001, y=.001, pmax=12, g1=4 (historic Dirac+anti convention). In particular **mφ=.001 is NOT mφ=86 eV**; no physical Yukawa parameter fitting is implied. The analytic Γ_rest=2.74969609752×10⁻⁹ in these arbitrary energy units. For f2=FD(T2=1) and bφ=BE(Tφ=1), different prescribed f1 temperatures give **opposite** physical net φ production:

| Tχ1 | p radial n8 Rφ | n16 Rφ | n32 Rφ | n48 Rφ | H at n48 |
|---:|---:|---:|---:|---:|---:|
| 0.9 | −2.656292×10⁻¹¹ | −2.658201×10⁻¹¹ | −2.658241×10⁻¹¹ | −2.658241×10⁻¹¹ | +7.419917×10⁻¹² |
| 1.0 | ≈0 | ≈0 | ≈0 | ≈0 | ≈0 |
| 1.1 | +2.736114×10⁻¹¹ | +2.738074×10⁻¹¹ | +2.738107×10⁻¹¹ | +2.738107×10⁻¹¹ | +6.934143×10⁻¹² |

**Negative for simplified “always produces φ” assumptions:** when χ1 is colder than the common χ2/φ bath at these synthetic conditions, the reverse χ2+φ→χ1 reaction dominates, **destroying** net φ; the quantum H theorem is positive in BOTH directions. This is not a new cosmic relic result, but proves that keeping inverse decay is essential for full evolution. The single Gauss angular n16 is enough for these smooth thermal shapes, not for arbitrary Bose IR spectral gradients.

Independent chemical fixed point test at T=.8, μ1=.03, μ2=.05, μφ=−.02 (nonzero and μφ<min Eφ) returns |Rφ|≈4.88×10⁻²⁷ and entropy≈1.21×10⁻⁴¹ at nrad24/nangle16. Eventwise energy defect is below 2.3×10⁻¹⁵ for tested moving mothers, total weak Q1+Q2+Qφ below ~4×10⁻²⁶ from floating roundoff. Four focused unit tests passed locally after fixing a genuine vectorized-angle validation bug; that original short-lived test failure was corrected before commit, not hidden by weakening physical criteria.

## Reproducibility, negative scope and next stage

[Executable moving Yukawa Pauli/Bose weak source](moving_yukawa_bose_pauli_weak.py), [4 tests](../tests/test_moving_yukawa_bose_pauli_weak.py), [machine data](results/moving-yukawa-bose-pauli-weak-2026-11-22.json), [permanent stage](../archive/stages/2026/11/2026-11-22-moving-yukawa-bose-pauli.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_moving_yukawa_bose_pauli_weak.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.moving_yukawa_bose_pauli_weak

This is an **independent Lorentz-boost / angular-integral extension** of earlier project rest-width and pointwise quantum controls, not rebranding the old paper. Method context: [Yin 2023](https://arxiv.org/abs/2301.08735) and [Sakurai & Yin 2025](https://arxiv.org/abs/2410.18968); no external source code was copied.

**Missing P0:** paired **mode-resolved** χ1,χ2,φ strong source with finite physical density-of-states cells (φ Bose!), no negative interpolation, exact three-species energy and χ number plus quantum H, and time evolution actually producing χ2 and real φ from dynamically depleted χ1 while simultaneously scattering in the ONE lepton bath. Need independent 3+ moving parent/decay daughter energy/angular/IR/cutoff convergence and actual χ/φ gauge/SM/Zprime+FRW. The 86eV benchmark, Γelastic/H, φ freezeout, absolute ΩDM, physical BEC and CLASS-ready remain **unvalidated**. Scientific status in_progress.
