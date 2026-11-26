# Rigorous near-rest 86-eV thermal Yukawa decay absolute rate upper bound (2026-11-24 UTC)

## Scientific question and restricted scope

A χ1→χ2+φ event producing an exactly rest-frame 86-eV scalar requires an extremely high-energy parent (p*=55,813.953488 GeV) in the illustrative tree-level Yukawa model. Can the actual *Bose-enhanced forward and inverse* thermal 1↔2 collision source in a finite near-rest energy window be rigorously bounded, without underflow or assuming Bose amplification is absent?

Yes, **only** for explicitly fixed Fermi-Dirac χ1/χ2 and strictly subcritical Bose-Einstein φ occupations at stated temperatures and chemical potentials. This is an analytically controlled rate **upper bound**, not a computed collision rate, nonlinear BEC no-go, FRW history or observational constraint.

Masses (GeV): M=.1, m2=.02, mφ=86e−9; Yukawa y=.001, g1=4, isotropic inclusive rest width Γrest. For T1=.9, T2=Tφ=1 GeV and μ1=μ2=μφ=0 we quantify four daughter windows.

## Exact parent support and mathematical upper bound

For any scalar laboratory k∈[mφ,K], with K<Eφ*≈.048 GeV, the exact daughter-parent rapidity theorem yields the finite **nested union** p∈[a,b], where

- ξφ=asinh(kstar/mφ), ζK=acosh(K/mφ);
- a=M sinh(ξφ−ζK), b=M sinh(ξφ+ζK).

This requires no finite parent UV cutoff: the on-shell parent momenta emitting a massive scalar *within this near-rest window* are intrinsically bounded both above and below.

The actual signed quantum event occupation difference is J=F−R,
F=f1(1−f2)(1+bφ), R=(1−f1)f2 bφ. Since 0≤f1,f2≤1 and bφ≥0,

    |J|≤F+R≤f1(1+bφ)+f2*bφ.

For FD populations with T>0, fFD(E;T,μ)≤exp[−(E−μ)/T]. Every allowed event has E1≥p≥a and E2=E1−k≥a−K. For BE at chemical potential μφ<mφ, bφ(k)≤bmax=1/expm1((mφ−μφ)/Tφ), and therefore

    |J|≤exp[−(a−μ1)/T1]*(1+bmax)
         +exp[−(a−K−μ2)/T2]*bmax.

The **same positive physical inclusive parent event measure** is
dW=g1 Γrest M p²dp/(2π² E1) * dc/2. The allowed angle subset has measure ≤1 and p²/E1≤p. Thus the integrated absolute signed scalar number source over the full near-rest daughter energy window satisfies

    |Rφ[mφ,K]| ≤ g1 Γrest M(b²−a²)/(4π²) *
        { exp[−(a−μ1)/T1]*(1+bmax)
         +exp[−(a−K−μ2)/T2]*bmax }.

The absolute scalar energy source obeys |Qφ[mφ,K]|≤K times the same right-hand side. These are strict conservative **upper bounds** within the stated trial thermal model and positive mass on-shell kinematics; no huge parent momentum integration is required. Log10 implementation uses stable two-term log-sum to avoid false underflow zeros.

## Numerical bounds in the prescribed thermal benchmark

Parameters: T1=.9,T2=Tφ=1 GeV, μ1=μ2=μφ=0, and thermal BE occupation at mφ bounded by bmax≈11,627,906.4767.

| Scalar window | Exact union of physical parent momenta (GeV) | log10(|Rφ| bound / GeV⁴) | log10(|Qφ| bound / GeV⁵) |
|---|---|---:|---:|
| 86–86.1 eV | [53,186.4833,58,571.2236] | **−23093.3062** | **−23100.3712** |
| 86–87 eV | [47,926.6978,64,999.2081] | −20808.5047 | −20815.5652 |
| 86–96 eV | [34,616.6849,89,991.2113] | −15027.4858 | −15034.5035 |
| 86–186 eV | [13,678.1666,227,749.6322] | −5933.1283 | −5939.8588 |

At 86–86.1 eV the **absolute rate bound** is ~10^(-23093.3062) GeV⁴. This is overwhelmingly smaller than the O(10^-11) GeV⁴ total weak rate of the parent 1-GeV-scale thermal benchmark, even after including the physical 1+bφ Bose factor. For μφ=mφ−1e−16 GeV, bmax>10^15 but the upper bound remains <10^-23000 GeV⁴. At the exactly critical μφ=mφ this proof's uniform BE cap is divergent and no finite bound is asserted.

**Scientific interpretation:** the *isolated on-shell 1↔2 channel* at these fixed modest temperatures cannot directly feed near-rest 86-eV scalar excited states at an appreciable thermal rate. This does not exclude nonthermal Bose amplification, ultrahigh temperatures, number-changing/2→2/SM channels, cosmological redshift, a singular delta zero mode, a true condensate, or a different microscopic theory. A finite upper bound on a continuous near-rest excited-state interval cannot be naively applied to a separate condensate δ³(p).

## Reproducibility, negative limits and next action

Source: [phi_bose_near_rest_thermal_bound.py](phi_bose_near_rest_thermal_bound.py).
Tests: [test_phi_bose_near_rest_thermal_bound.py](../tests/test_phi_bose_near_rest_thermal_bound.py).

Run commands:
- OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_near_rest_thermal_bound.py
- python -m kinetics.phi_bose_near_rest_thermal_bound

Four focused local tests passed in 0.47s: exact Lorentz parent interval nesting/direct boost support; log10 absolute bound regression; monotone widening and near-critical thermal Bose occupation; exclusion of critical μφ from this theorem. New public CI separately verify.

**NEXT_ACTION:** register this bound as a *conditional physical infrared gate* in the coupled kinetic plan. Avoid attempting direct brute-force GeV-temperature thermal parent quadrature at ~55TeV unless the studied occupations/processes genuinely change. Continue development of a **fully coupled physical finite-time χ1/χ2/φ strong collision method** with common positive event weights, real DOS capacities, two charges/energy, exact quantum detailed balance, Pauli/Bose bounds and H theorem. No SM+FRW/86-eV relic/CLASS or galactic claims.
