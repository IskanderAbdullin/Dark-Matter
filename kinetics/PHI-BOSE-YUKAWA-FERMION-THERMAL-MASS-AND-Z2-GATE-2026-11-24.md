# Physically decisive thermal-mass gate for the 86-eV real Yukawa scalar: exact two-Dirac static one-loop potential (2026-11-24 UTC)

## Question, assumptions, and why the old vacuum-mass kinetic extrapolation is incomplete

The existing three-field χ1↔χ2+φ physical on-shell code uses **fixed vacuum masses** m1=.1 GeV, m2=.02 GeV, mφ=86e−9 GeV and off-diagonal Dirac Yukawa y=.001. Earlier physical Bose IR accessibility estimates (including a >55-TeV parent momentum to hit rest φ) use this **VACUUM** dispersion. Before embedding them in a plasma FRW history we must independently check a minimal radiative thermal correction to the φ dispersion and to its candidate stabilizing Z2.

Take a *specific conditional* microphysical model: two **Dirac** fermions χ1,χ2 at common thermal Fermi-Dirac temperature T, **zero chemical potentials**; a real scalar φ with off-diagonal real Yukawa vertex

    L_int=-y φ (barχ2 χ1+barχ1 χ2),
    M_F(φ) = [[m1,yφ],[yφ,m2]],   m1=.1GeV,m2=.02GeV.

A candidate Z2 can send φ,χ1→minus themselves, leaving χ2 unchanged. Thus V_F(T,φ) is even in φ. Only a **one-loop free equilibrium fermion thermal effective-potential curvature about φ=0** is being evaluated. It is **not** automatically the real-time scalar pole/dispersion or a prediction of a cosmological phase transition: actual occupation and thermal/SM self-energies, daisy resummation, φ self-interaction and the finite dark+SM history are not supplied.

## Independent derivation with all spin multiplicities

For each Dirac eigenmode with field-dependent mass eigenvalue m_i(φ), its thermal free-energy potential at μ=0 is

    V_F(T,φ) = -4 T ∫ [d³p/(2π)³] Σ_{i=1,2} ln[1+exp(-sqrt(p²+m_i(φ)²)/T)].

The factor four counts two spin and particle/antiparticle states of EACH Dirac fermion. At φ=0,

    dV_i/d(m_i²) = 2 I_F(m_i,T),
    I_F(m,T)=∫ d³p/(2π)³ [1/(E_p(exp(E_p/T)+1))].

Since Tr[M_F(φ)²] = m1²+m2²+2y²φ² exactly, the high-T mass correction is positive and equals **δmφ,F²→ y² T²/3** for two thermally populated Dirac modes. At general finite m/T, the off-diagonal mass matrix yields

    d²m_1²/dφ²|0=+4m1y²/(m1−m2),
    d²m_2²/dφ²|0=−4m2y²/(m1−m2),

hence the EXACT one-loop finite-T static curvature

    δmφ,F²(T)= [8 y²/(m1−m2)] [m1 I_F(m1,T) − m2 I_F(m2,T)].

This expression has **no free BGK rate, no arbitrary pole-mass ansatz**, and is computed by a one-dimensional positive Fermi integral for each Dirac eigenstate. Its sign is not fixed at intermediate T: with m2≪m1 the lighter fermion can dominate the negative mass-mixing response while the heavier is Boltzmann-suppressed.

This general mass-dependent expression was independently checked in the Python tests against a **separate, direct 2×2 eigenmass determinant integration V_F(φ,T)** and its central second derivative (two background fields ±φ), without reusing the analytical I_F expression. At T=.01,.02,.035,.04,.1,1GeV the analytic/finite difference relative discrepancy is ≤8e−4 with φ=.2GeV, dropping with φ² as expected.

The high-T limit is also checked against the universal finite-T Dirac fermion potential expansion; references: [Dolan & Jackiw (1974)](https://doi.org/10.1103/PhysRevD.9.3320), [Quirós, Finite Temperature Field Theory and Phase Transitions](https://arxiv.org/abs/hep-ph/9901312). These papers provide the underlying thermal field-theory framework, **not** an independent measurement or validation of this particular off-diagonal model.

## Numerically converged conditional result: negative *static curvature* window

Let the vacuum φ mass-squared be mφ,0²=(86e−9GeV)²=7.396e−15GeV². Without additional scalar quartic or portal, the quadratic curvature of V_eff at φ=0 is mφ,0²+δmφ,F²(T). Finite-mass Fermi quadrature gives:

| T (GeV) | δmφ,F² (GeV²) | Total static curvature (GeV²) | Meaning |
|---:|---:|---:|---|
| .002 | −7.557944e−17 | +7.320421e−15 | locally positive |
| .005 | −1.257047e−13 | **−1.183087e−13** | locally negative |
| .010 | −2.706703e−12 | **−2.699307e−12** | negative |
| .020 | −1.79132110e−11 | **−1.79058150e−11** | negative |
| .030 | −1.923808e−11 | **−1.923068e−11** | negative |
| .035 | −2.439809e−12 | **−2.432413e−12** | negative |
| .040 | +3.044558e−11 | +3.045298e−11 | locally positive |
| .100 | +1.944042e−9 | +1.944050e−9 | positive |
| 1.000 | +3.291148e−7 | +3.291148e−7 | positive |

Two roots of **the strictly specified static curvature** occur at T=**0.00324536695GeV** and **0.03548213975GeV** (≈3.245 and 35.482 MeV). Between them the φ=0 quadratic potential is negative **if** both Dirac modes really have equilibrium thermal occupations at the same T and scalar self interactions/portals are absent. This raises a **conditional Z2-stability/reentrant-symmetry-breaking concern**, not an established phase transition. Real thermodynamic minima, scalar self-interactions, collision freeze-out, finite-time nucleation, vacuum renormalization and in-medium masses have NOT been computed. The vacuum Z2 is not thereby proven broken in the real Universe.

## Scalar self-coupling is a critical missing physical parameter

Even a small bounded scalar quartic V_self=λφ⁴/4! contributes leading **positive** static thermal self-energy λ T²/24 for a light real scalar. In this simple perturbative one-loop approximation, requiring total curvature nonnegative **at every T in the examined 0.004–0.035 GeV interval** yields

    λ >= max_T [-24(mφ,0²+δmφ,F²(T))/T²]
      ≈ **1.08515714×10⁻⁶**,
    at T≈**0.0185984949 GeV**.

This is only a **conditional sufficient static-curvature criterion** in the specified truncated potential at leading order, not a global vacuum/phase-history proof; a quartic can also influence φ self-scattering, number-changing rates, and perturbation dynamics, which must be included before assuming a stable 86-eV cold condensate. The no-quartic numerical prototype therefore omits an important parameter which can affect even the local potential stability.

## High-temperature static thermal scale: not the original 86eV dispersion

At T=1GeV, exact finite-mass integral gives δmφ,F²=**3.29114832e−7 GeV²**, corresponding to a **STATIC curvature mass scale** sqrt(mφ,0²+δm²)=**0.000573685 GeV=574keV**, a factor ≈6670 larger than the vacuum 86eV. At T=100GeV, the high-T coefficient is reached to ~3e−6 relative and sqrt curvature ≈0.05773GeV. Replacing only the scalar mass by mφ,eff(T)≈yT/sqrt3 in vacuum 1↔2 kinematics would make mφ,eff≈m1−m2=.08GeV near T≈**138.564GeV**, but **this phi-only insertion is not a consistent claim that the actual thermal 1↔2 channel closes there**: the fermion dispersions, spectral cuts, thermal widths, SM bath, in-medium mixings and full quantum self-energy all need to be computed together. Static effective-potential curvature and real-time on-shell pole mass are distinct observables.

This is a stronger, model-specific reason that the existing vacuum φ-rest parent threshold and fixed-GeV narrow Bose IR source conclusions cannot be extrapolated unchanged to early thermal epochs. An 86-eV scalar at T=0 is not necessarily an 86-eV quasiparticle in a populated plasma.

## Reproduction, negative evidence and next stage

Code: [two-Dirac thermal I_F and direct determinant](phi_bose_fermion_thermal_mass.py); [five independent tests](../tests/test_phi_bose_fermion_thermal_mass.py); [numerical table](results/phi-bose-fermion-thermal-mass-2026-11-24.json); [dated permanent stage](../archive/stages/2026/11/2026-11-24-yukawa-fermion-thermal-scalar-curvature.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_fermion_thermal_mass.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_fermion_thermal_mass

**NEXT_ACTION:** derive finite-T *real-time* on-shell scalar and χ1/χ2 self-energies/dispersion (including HTL, finite-momentum dependence), select and justify λφ, thermal scalar Z2 history and SM portal/chemical equilibration conditions. Update 1↔2 and 2↔2 rates with the SAME in-medium kinematics, physical degeneracies, energy/H entropy and FRW source; then obtain a real frozen fφ(q), Ω and observational transfer. Do not declare a physical thermal phase transition, stable BEC, or cosmological particle production from this static curvature alone.
