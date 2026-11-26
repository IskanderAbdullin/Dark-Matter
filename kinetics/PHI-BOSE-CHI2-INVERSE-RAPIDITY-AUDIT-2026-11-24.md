# Massive fermion daughter χ2: **two-sided** inverse parent support and third independent Yukawa N/Q integration order (2026-11-24)

## Scientific issue

Previous physical Yukawa 1↔2 weak and conditional scalar-mode strong checks integrated parent momentum first, or scalar φ energy k first. For a fixed **massive** daughter χ2 energy ε2, however, one generally cannot replace physical on-shell parent integration by [a minimum parent momentum, P]. The allowed parent momentum can have a **finite upper bound strictly below P**, and naive parent Gauss quadrature can badly miss the Bose-enhanced infrared. This is a kinematic effect, not a condensate or physical freeze-out.

Benchmark m1=.1, m2=.02, mφ=86e−9 GeV, Yukawa y=.001, inclusive isotropic unpolarized Γrest, g1=g2=4 illustrative internal occupation counts. The channel is χ1↔χ2+φ, actual massive daughter dispersion and physical p²dp density of states. Parents p∈[0,P], P=2/8/16 GeV.

## Exact relativistic inverse rapidity theorem

Define rest-daughter rapidity ξ=asinh(kstar/m2), lab daughter rapidity ζ=asinh(sqrt(ε2²−m2²)/m2), and parent rapidity y=asinh(p/m1). The Lorentz boost yields ε2=m2 cosh ξ cosh y ± m2 sinh ξ sinh y in the extreme collinear directions. Its **full** on-shell energy shell condition is

    |ξ−y| <= ζ <= ξ+y

which is equivalent to the *two-sided* allowed parent rapidity

    |ξ−ζ| <= y <= ξ+ζ.

After hard parent pmax=P, the true momentum support is

    p_low=m1 sinh|ξ−ζ|,
    p_high=m1 sinh[min(asinh(P/m1),ξ+ζ)].

Return no support if the lower rapidity exceeds asinh(P/m1). This proof is independent of the Yukawa matrix-element occupation factors; rates using uniform rest-angle additionally assume isotropic unpolarized decay.

At daughter χ2 at rest, ε2=m2=.02, parent p is a single special value p*=m1*kstar/m2≈**0.24 GeV**. At ε2=.0201 GeV and P=2, the actual interval is **[.2151675203,.2672324797] GeV**. At the parent-rest-frame daughter energy ε2*=.052 GeV, interval approximately **[0,1.248] GeV**, not [0,2]. An algorithm integrating every admissible daughter mode up to P would add forbidden events and corrupt weak normalization.

## Exact physical daughter-mode generator and an efficient variable change

The physical daughter DOS per energy is g2 ε2 sqrt(ε2²−m2²)/(2π²). Dividing the same positive on-shell event measure by this DOS gives, with χ1 and φ held fixed,

    dot f2(ε2) = G2(ε2)[1−f2(ε2)] − L2(ε2) f2(ε2),

where G2 is the positive parent integral of f1(E1)[1+bφ(E1−ε2)] and L2 the positive parent integral of [1−f1(E1)]bφ(E1−ε2); both use prefactor g1 Γrest m1²/[2g2 kstar ε2 sqrt(ε2²−m2²)] times p dp/E1. This guarantees exact conditional Pauli bounds 0<=f2(t)<=1 under frozen other distributions with the exponential propagator f2(t)=G/(G+L)+(f2(0)−G/(G+L))exp[−(G+L)t]. Common-T chemical equilibrium recovers the exact FD occupancy at μ2=.05, μ1=−.15, μφ=−.2 (tested ε2=.03/.052/.1/.4 GeV).

Crucial conditioning: **at fixed ε2, k=E1−ε2, hence p dp/E1=dk**. Therefore integrate the two-sided parent interval using a positive *logarithmic scalar-energy k* Gauss map, rather than a broad fixed p Gauss. Then integrate the outer χ2 rapidity ζ (physical χ2 DOS and dε2=p2 dζ), splitting at changes of the support interval. This is independent of the previous parent-first and φ-energy-first implementations.

The same physical event weights, multiplied by actual parent energy E1, give **Rφ, Qφ, Q1, Q2** with Qφ=(parent-energy event −χ2 energy event); exact event-level E1=ε2+k is respected.

## Quantitative three-way check, negative coarse-grid evidence

Max relative error over all four moments Rφ,Qφ,Q1,Q2 against the independent analytic log/dilog thermal reference:

| Parent P GeV | daughter ζ and inner log-k 8×8 | 64×64 | adaptive outer ζ, inner log-k 36 |
|---:|---:|---:|---:|
| 2 | **9.28e−3** | 2.53e−5 | **≤4.5e−13** |
| 8 | 2.53e−4 | 4.02e−6 | **≤3.6e−13** |
| 16 | 6.55e−3 | 8.80e−9 | **≤3.1e−10** |

The first attempted **untransformed** p-Gauss and daughter ζ Gauss at P16 gave ~53.3% four-moment relative error at 8×8, still ~7.8% at 64×64: a serious negative quadrature finding. Exact p/E1 dp=dk plus positive log-k weights greatly reduces that error; adaptive outer rapidity resolves residual physical support and IR cutoffs. No wrong-result tolerance has been hidden.

At P2 the new independent χ2-rapidity adaptive result Rφ=−1.29545395549e−11 GeV⁴ matches both the independently validated φ cut-cell and the analytic 1D parent references. All three energy ledgers agree within adaptive error; independently integrating each Q introduces a small P16 numerical summed Q defect ~1.23e−20 GeV⁵, **not** an exact numerical conservation claim. A future coupled strong method must use one shared event extent for strict discrete energy conservation.

Four focused local pytest tests passed in 4.28 seconds, checking finite two-sided support against independently boosted on-shell energy inequalities, exact Pauli semigroup and μ equilibrium, independently convergent 4 moments at P2/8/16, and explicit coarse grid IR misses.

## Limits and NEXT_ACTION

Source: [physical χ2 inverse rapidity integration](phi_bose_daughter_inverse_rapidity.py). Tests: [four focused regressions](../tests/test_phi_bose_daughter_inverse_rapidity.py). Reproduce with OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_daughter_inverse_rapidity.py and OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_daughter_inverse_rapidity. Public Scientific CI separately verify after commit.

**Not a solved coupled three-species strong solver:** χ1, φ held fixed during χ2 conditional mode integration; no continuous time energy conservation with all species evolving, no finite SM/FRW, true Bose zero mode, Ω or CLASS. Daughter-state degeneracy g2 must be matched to an explicit microscopic spin/channel model before a complete physical prediction.

**NEXT_ACTION:** use the same general two-sided inverse parent support for scalar φ as well (the φ-at-rest required parent may be enormous), enforce physical daughter DOS and Pauli/Bose mode bounds in a fully coupled common-event nonlinear operator, prove both dark charges and all three Q/H, then independently converge pcut/qmin/dt/true thermal bath.
