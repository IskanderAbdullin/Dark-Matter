# Independent massive-lepton 5D χe heat benchmark — 2026-11-22

## Hypothesis and pre-registered failure gates

Test the **conditional broken-electroweak** Dirac-vector contact interaction
G(χ̄γμχ)(ēγμe), with distinct massive fermions, full four-Fermi
gain/loss and zero chemical potentials. Does the proposed **direct**
one-dimensional momentum-transfer kernel reproduce the *independently
constructed* five-dimensional center-of-mass (CM) scattering heat integral,
including finite electron mass? Failure criteria: a persistent relative
discrepancy >1e-3 after angular/radial/transfer refinement, nonzero common-T
heat beyond rounding, wrong sign of heat for Tχ<Te versus Tχ>Te, or a
Dirac-spin trace mismatch. Numerical agreement is **not** a physical
γ_el, an eigenvalue, or a 86-eV dark-matter prediction.

Independent checks: (i) explicit 4×4 Dirac matrices; (ii) tensor-product
Gauss quadrature over p,k,μ,cos θ*,φ* with on-shell CM→laboratory boost,
**without invoking the q kernel**; (iii) separate E,ω,Δ=E′−E>0,ℓ
energy-transfer quadrature using the claimed q kernel; (iv) massless
5D independent reference obtained earlier.

## Exact formula tested

For M=mχ, m=me, s=M²+m²+2(Eω−pkμ),
λ(s,M²,m²)=[s−(M+m)²][s−(M−m)²],
p*²=λ/(4s), t=−2p*²(1−cos θ*), D=s−M²−m²,
u−M²−m²=−D−t. Explicit gamma-matrix traces give

    |M|²_avg = 2G²{ D² + (D+t)² + 2t(M²+m²) }.

With F=fχ(E)fe(ω)[1−fχ(E′)][1−fe(ω′)] and
R=fχ(E′)fe(ω′)[1−fχ(E)][1−fe(ω)], the independent COM heat
is the **full-domain** integral

    Qχ = gχ ge/(512π⁵) ∫dp p² ∫dk k² ∫dμ ∫dc ∫dφ/(2π)
         × vM / s × |M(s,t)|² × (E′−E) × [F−R],
    vM = sqrt[(Eω−pkμ)²−M²m²]/(Eω).

The 1/2 relative to the forward-only collision measure is required here
because both time-reverse orientations are integrated. The outgoing χ
energy is obtained from a genuine CM boost:
E′=γ(Eχ*+βp* cosα′), Eχ*=(s+M²−m²)/(2√s);
cosα′=cosα cosθ*+sinα sinθ* cosφ,
cosα=γ(p_parallel−βE)/p*, with
p_parallel=p(p+kμ)/|p+k|. This approach does not use momentum-transfer
cones or an angle-averaged trace.

Separately, the energy-transfer **half-domain** integral is

    Qχ = ∫_M^∞dE ∫_m^∞dω ∫_0^(ω−m)dΔ
         Δ [F−R] K(E,ω;E+Δ,ω−Δ),

where K=gχ ge/(256π⁵) ∫_(ℓmin)^(ℓmax) dℓ〈|M|²〉cone
and ℓmin=max(|p−p′|,|k−k′|),
ℓmax=min(p+p′,k+k′). There is **no** additional 1/2 because
only Δ>0 is counted. Under common T and the same μχ,μe on both
sides, R=F exactly; for Tχ<Te the heat must flow into χ.

The independent 4×4 gamma trace agrees with the analytic massive-spin
polynomial at four physical (M,m,s,cosθ*) points to <=1.5e-15 relative
in the local check (see executable
tests/test_finite_mass_dirac_trace.py). This is an algebraic check,
distinct from the phase-space normalization check.

## Fully specified finite-grid experiment

To avoid extraordinary physical-rate scales, numerical integration
used m=1 as an **energy unit**, M/m=0.1, G m²=0.001,
gχ=ge=4, Tχ/Te=0.9 or 1.1, and Te/me=0.1,1,10.
Both species were assigned thermal Fermi-Dirac occupation with μ=0.
E_in,max=M+30 Te and ω_in,max=m+30 Te; outgoing legs were not
independently hard-clipped, so the exponentially small cutoff tail
is an additional systematic. The normalized Q/(G² me⁹) is invariant
under physical mass rescaling. One direct rescaling check with
me=0.51099895 MeV and G=1 GeV^-2 at Te/me=1 yielded
Q_5D≈3.97897663e−32 GeV⁵ (n32), matching the dimensionless
normalization below. These are **not physical early-universe values**:
the contact UV completion and electron/positron evolution are missing.

| Te/me | Tχ/Te | independent 5D Q/(G²me⁹), n40 | direct ℓ n160 / 5D n40 − 1 |
|---|---:|---:|---:|
| 0.1 | 0.9 | +1.50804492085632e−14 | +5.69e−7 |
| 0.1 | 1.1 | −4.13146304202007e−14 | +8.11e−7 |
| 1 | 0.9 | +1.67490950359610e−2 | −7.05e−7 |
| 1 | 1.1 | −3.78116454982979e−2 | −5.32e−7 |
| 10 | 0.9 | +1.98565557482877e7 | +7.10e−10 |
| 10 | 1.1 | −4.43192337435731e7 | −6.61e−7 |

Full 5D uses tensor-product Gauss orders 24,32,40 per incoming radial
and angular dimension; transfer uses n=64,96,160 for each of E,ω,Δ
and nq=24 for ℓ. n40 and n160 show <1e−6 agreement in all six
**specific** benchmarks. The earlier massless independent 5D reference
Qχ=+1.9851492233308662e−8 GeV⁵ for M=.1GeV, me=0, Te=1GeV,
Tχ=.9GeV, G=.001GeV−², cut22GeV was reproduced at order 18
within 2.08e−6 relative.

**Negative precision result:** refinement is NOT monotone. At Te/me=1
and Tχ/Te=.9, transfer n64 gives +1.17e−5 relative to n40 5D,
n96 gives −2.46e−6, n160 −7.05e−7. At fixed n96, the cutoff
25→30→35 Te shifts Q by about 9e−6 relative, substantially
larger than the ℓ quadrature nq24→36 shift (~4e−9 relative).
Thus a coincidentally tiny n160 residual **does not certify**
uniform 1e−6 global error or physical convergence. A sharp
kinematic support boundary and exponentially weighted high-energy
tails remain to be bounded adaptively.

Data: [finite-mass-5d-independent-2026-11-22.json](results/finite-mass-5d-independent-2026-11-22.json).
Reproduction (each command is one independent case):

    python -m kinetics.independent_massive_chie_5d --temperature-ratio 0.1 --chi-temperature-ratio 0.9
    python -m kinetics.independent_massive_chie_5d --temperature-ratio 1.0 --chi-temperature-ratio 1.1
    python -m kinetics.independent_massive_chie_5d --temperature-ratio 10.0 --chi-temperature-ratio 0.9
    python -m pytest -q tests/test_independent_massive_chie_5d.py tests/test_finite_mass_dirac_trace.py

The other three combinations use the analogous temperature ratios.
Run --quick for CI/smoke. Tests cover the independent prior massless
reference, massive heat sign and numerical agreement, common-T
equilibrium and an independent explicit Dirac trace. No arbitrary
posthoc W averaging or BGK rate is involved.

## Literature and scope

Covariant fermionic collision operators and their Pauli treatment:
C. Pitrou, *Kinetic theory of fermions in curved spacetime*
(https://arxiv.org/abs/1701.08844).
For the independent missing electron-positron-chemistry gate,
*Electron-positron annihilation freeze-out in the early universe*,
Phys. Rev. D 101, 063507 (2020), DOI:10.1103/PhysRevD.101.063507.
These sources motivate the collision and plasma ledger, **not**
a validation of the present hypothetical χ portal.

**Still unproven:** UV/gauge-invariant mediator and finite-T masses,
real e± annihilation and neutrino/photon finite heat, continuum
nonthermal relaxation spectrum and Pauli-safe coupled evolution,
thermal χ1↔χ2+φ Bose production, H/gstar, frozen f(q), macroscopic
zero mode/BEC, absolute ΩDM, 86-eV physical relic, CLASS/P(k).
PR #1 must remain draft; CLASS-ready=false.
