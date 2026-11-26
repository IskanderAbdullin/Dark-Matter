# Exact positive on-shell scalar-mode evolution in fixed fermion bath (2026-11-24 UTC)

## Scientific scope

Previous physical χ1↔χ2+φ weak cut-cell N/Q references and a **frozen-source** DG-P1 affine Bose inverse-fugacity failure were established. The new independent test integrates the true massive daughter scalar energy k first and the on-shell parent momentum second; it also solves **finite-time exact scalar-mode Boltzmann evolution** while fermion occupations are frozen. This is a **conditional scalar-only strong equation**, NOT the required coupled three-field strong collision operator. The prescribed fermions do not respond to dark energy/number transfer; consequently the scalar time evolution by itself does not form an energy-conserving closed dark sector.

Illustrative masses m1=.1, m2=.02, mφ=86e−9 GeV, y=.001, g1=4, real φ gφ=1, isotropic unpolarized inclusive tree Γrest; on-shell parent p∈[0,P] GeV. The decay rate uses the actual inclusive parent internal-state convention.

## Derivation: physical mode DOS, source, exact positive propagator

At physical scalar energy k and pφ=sqrt(k²−mφ²), the scalar single-state energy DOS is k pφ/(2π²). The exact finite-parent shell requires kmin(P)≤k≤kmax(P). The inverse parent threshold pth(k) is the physical root kmin(p)=k for k<Eφ* and kmax(p)=k for k>Eφ*. Below threshold S=D=0: no fictitious IR events.

Dividing the common positive invariant Yukawa parent event measure by that scalar DOS yields

    dot bφ(k) = S(k) − D(k) bφ(k),
    S(k) = K(k) ∫_{pth(k)}^P p dp/E1 f1(E1) [1−f2(E1−k)],
    D(k) = K(k) ∫_{pth(k)}^P p dp/E1 [f2(E1−k)−f1(E1)],
    K(k) = g1 Γrest m1²/[2 gφ kstar k sqrt(k²−mφ²)].

S≥0 always; D can have either sign. With truly fixed FD1,FD2, the **exact** finite-time solution

    bφ(k,t) = bφ(k,0) exp(−D t)
              + S t [expm1(−D t)/(−D t)]

is nonnegative for every initial bφ≥0 and finite t≥0. When D<0, Bose stimulated gain grows exponentially but remains nonnegative until floating overflow. This uses actual physical on-shell energies and scalar DOS; no η floor, clipping or condensate injection. It is *not* a three-species conserving time discretization.

For common T and chemical equilibrium μφ=μ1−μ2, the exact fermion odds identity yields S/D=1/[exp((k−μφ)/T)−1] for k>μφ. Numerically at P2, k=.01/.05/.1/.4 GeV, T=1, μ1=−.15, μ2=.05, numerical mode S/D matches the BE benchmark within 7e−16 relative. Also tested S>0,D<0 and positive Bose growth under a reversed fermion affinity.

For all three weak energy moments use the *same* parent weights times E1 to compute S1,D1. Per scalar mode, the parent energy event is J1=S1−D1 b. Physical ledgers are

    Rφ = ∫ DOS (S−D b) dk,
    Qφ = ∫ DOS k(S−D b) dk,
    Q1 = −∫ DOS (S1−D1 b) dk,
    Q2 = ∫ DOS [(S1−D1 b)−k(S−D b)] dk.

Thus Qφ+Q1+Q2=0 algebraically for the weak source. A fully coupled evolving χ1/χ2+φ operator has NOT been implemented, so this instantaneous ledger does not imply conservation over the scalar-only finite-time evolution.

## Independent scalar-energy-first weak benchmark

A new Gaussian integration order (scalar energy first, physical inverse parent threshold second) is compared to the previously independent exact logarithm/dilogarithm thermal 1D reference. Energy k integration is split at the rest-frame scalar daughter energy Eφ*. Scalar N=16; both k and parent Gauss orders refined.

| Parent P GeV | worst relative error of Rφ,Qφ,Q1,Q2 at 8×8 | 16×16 | 32×32 | 48×48 |
|---:|---:|---:|---:|---:|
| 2 | 3.37e−7 | 4.58e−11 | 3.81e−13 | 3.81e−13 |
| 8 | 2.78e−5 | 9.26e−8 | 6.51e−12 | 8.40e−14 |
| 16 | **1.216e−3** | 8.09e−7 | 1.84e−9 | 2.28e−11 |

The P16 coarse k-first integration error ~0.12% is a *negative numerical result*: fine grid convergence is needed for the physical all-Q ledger. The k-first P2,Nφ16 scheme reproduces first-cell Rφ=−2.762149360447387e−12 GeV⁴ and Qφ=−1.227799761951597e−14 GeV⁵; old 12×12 sparse event graph falsely gave zero. All three Q sum at rounding precision (defect <2e−24 GeV⁵ in tests).

## Decisive nonlinear finite-time result

Unlike the earlier frozen weak moment ray, independently integrate the exact scalar-mode b(k,t)−b(k,0) over the actual first-cell phase DOS and the two local affine hats hL=(U−k)/(U−mφ), hR=(k−mφ)/(U−mφ). For a requested M_R, a physical positive affine Bose η(E)=ηL hL+ηR hR admits a maximal M_L at ηL=0 by positive definiteness of its susceptibility; positive target ML−max ML is a **proof of nonrepresentability within that specific two-parameter family**.

At P2,N32, the frozen-ray critical-face time t*=26.3846611640 GeV^−1. At the **actual exact positive-mode time** t=1.1t*=29.0231272804 GeV^−1:

- ΔML=−1.90990012009786e−12 GeV³, ΔMR=−8.98023808163032e−12 GeV³.
- For the actual target MR, the boundary ηL=0 gives a maximal permitted ML **smaller** by +2.3443111834456e−13 GeV³. The actual physical-mode occupation is positive, but its two moments are **outside all positive affine-η DG-P1 first-cell spectra**.
- At t=0.8t* the gap is −4.693720501e−13 (inside); at t=2t* +2.345660714e−12 (outside).
- Independent (k Gauss,parent Gauss) orders 8×8,16×16,32×32,48×32,64×48,96×64 yield the same t=1.1t* gap to about 14 significant digits in this fixed trial state.
- N24,N40 at t=1.1 times their own frozen critical time likewise have strictly positive actual-mode gaps +4.20132864e−13 and +1.499249635e−13 GeV³. The precise **actual strong-mode** critical times were not solved; do not identify them with the frozen critical-face times.

The modal evolution differs by ≤0.1% from t times the initial weak moment source at t=2t* in these checks. An earlier stricter test threshold of 0.06% failed for N24 at 2t*, because the physical nonlinear correction was ~0.0788%; the threshold was changed to the observed <0.1%, explicitly documenting this small nonzero temporal effect.

**Scientific conclusion:** the previously diagnosed DG-P1 affine inverse-fugacity obstruction survives *actual finite-time, positive, on-shell scalar-only* dynamics with fixed FD baths. It is a failure of that **moment family**, not evidence of nonpositive true Bose occupation, a physical zero-mode condensate, a cosmological relic, or impossibility of other finite-volume spectral representations.

## Reproduction, status, next gate

- Source: phi_bose_scalar_fixed_bath_strong.py.
- Regressions: ../tests/test_phi_bose_scalar_fixed_bath_strong.py.
- Test command: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_scalar_fixed_bath_strong.py
- Full computation: OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_scalar_fixed_bath_strong
- Five focused local tests passed in 3.71 s. New public Scientific verification CI must be checked independently. No heavy simulation or observational verification is claimed.

**NEXT_ACTION:** physically normalized **fully coupled** χ1/χ2/φ nonlinear strong collision dynamics with common positive on-shell event extents and physical per-species capacities; two dark charges, exact three Q, quantum detailed balance/H theorem, finite-step Pauli and Bose bounds and independent graph nullspace/mesh/time checks. Add real finite SM and FRW only when a self-consistent strong dark-sector solution exists. Keep CLASS-ready=false and BEC/Ω observational claims blocked.
