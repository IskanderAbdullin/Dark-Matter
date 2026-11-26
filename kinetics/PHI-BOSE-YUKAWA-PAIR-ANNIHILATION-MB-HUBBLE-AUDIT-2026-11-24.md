# Off-diagonal Yukawa fermion-pair → two real scalars: vacuum matrix element and conditional forward MB rate (2026-11-24)

## Motivation: are the missing 2↔2 reactions physically negligible?

The existing physical-capacity quantum χ1↔χ2+φ finite-time solver is **closed under one 1↔2 channel only** and therefore conserves two positive number combinations N1+N2 and N1+Nφ. The [exact full-network stoichiometric rank audit](PHI-BOSE-YUKAWA-CONTINUOUS-CHARGE-VS-ACCIDENTAL-COLLISION-2026-11-24.md) showed that this conservation need not extend to other reactions allowed by the *same* renormalizable off-diagonal Yukawa vertex. To decide whether those pair channels might matter compared with Hubble expansion, compute an actual **spin-averaged vacuum tree cross section** and independently normalized **forward-only Maxwell-Boltzmann reaction density** for

    χ2 + antiχ2 -> φ + φ     (virtual χ1 in t/u channels)

with Dirac masses m1=.1GeV, m2=.02GeV, real scalar mφ=86e−9GeV, Yukawa y=.001. The final identical scalars carry 1/2! symmetry factor. The assumed toy Lagrangian is L_int=−y φ barχ2 χ1+h.c.; spin degeneracy gχ2=gantiχ2=2 for the **separate** particle and antiparticle populations. A possible Z2 odd φ,χ1/even χ2 permits this channel. This is a **model-conditional tree/MB calculation**, NOT a proven full χ+SM plasma kinetic rate, chemical relaxation, relic or BEC.

## Actual matrix element and an independently checked spin trace

For incoming p1(χ2),p2(antiχ2), outgoing k1,k2 identical φ, define qt=p1−k1,qu=p1−k2 and Dt=t−m1²,Du=u−m1². The real-scalar identical final-state t and u exchange **sum**:

    M=y² bar(vχ2)( (slash qt+m1)/Dt + (slash qu+m1)/Du ) uχ2.

The unpolarized initial spin average is

    |M|²_avg = y⁴/4 Tr[
       (slash p2−m2) A (slash p1+m2) gamma0 A† gamma0
    ],
    A=(slash qt+m1)/Dt+(slash qu+m1)/Du.

We implemented this in TWO genuinely independent ways: (i) a scalar invariant Minkowski-dot expansion of all four traces and (ii) direct numerical **4×4 complex Dirac gamma matrices and Hermitian Dirac adjoint** in the tests. For √s=.041,.05,.1,.4,1,2GeV and cosθ=0,.3,.9 the spin-averaged |M|² values agree to double precision. The amplitude is nonzero above threshold; near the nonrelativistic pair threshold χ2χbar2→φφ exhibits the expected p-wave behaviour σv_Møller∝β², not an unsuppressed s-wave.

With identical final φ bosons, the physical vacuum angle-integrated cross section is

    dσ/dcosθ = [1/(64πs)] * [pout/pin] |M|²_avg,
    σ(s) = ∫_−1^1 dcosθ dσ/dcosθ.

At √s=.05/.1/.5/1GeV the converged cross sections are respectively **1.44620002447e−12, 1.272097604216e−12, 8.58426302470e−14, 2.93535079494e−14 GeV^-2** for y=.001. This *already includes* initial spin average and final identical 1/2.

### Serious negative numerical control: naïve angle-Gauss can be wrong by ~46%

Physical t/u fermion propagators produce narrow *forward/backward* angular peaks when √s≫m1. A **36-node plain Gauss-Legendre cosθ** fails to resolve those peaks at high energy:

- √s=5GeV: coarse raw cosθ Gauss undercounts σ by **4.7993%**;
- √s=100GeV: it undercounts σ by **46.2870%**.

Initial MB rate estimates based on plain cosθ quadrature were therefore not reliable at high T. We fixed it by a **positive conformal angular map**

    cosθ=tanh u,    dcosθ=sech²u du,

with u integrated from −14 to14, Gaussian weights and 96/128-point refinements. This clusters nodes exponentially at each physical forward/backward peak, preserves positivity, and at √s=1/2/5GeV matches independent adaptive cosθ integration with **explicit very small absolute** error tolerance to roughly ≤1e−8 relative by 128 nodes. At √s=100GeV a mapped 96→128 refinement differs by <1e−5. The old 36-node ~46% error is independently recreated in the regression test, not simply excluded from reporting.

This error matters particularly for relativistic early-Universe parents. The present results use the **corrected** angular map.

## Forward Maxwell-Boltzmann rate versus Hubble: explicitly conditional

For *ideal, unblocked MB* χ2/antiχ2 populations at a specified T and with φ initially empty, take

    nχ2^MB = g2 m2² T K2(m2/T)/(2π²), g2=2,
    γ_forward = g2² T/(32π⁴)
                ∫_(4m2²)^∞ ds σ(s)(s−4m2²) √s K1(√s/T),
    Γ_forward/χ2 = γ_forward/nχ2^MB.

The integral is the invariant standard **MB reaction density** for two distinct fermion/antifermion initial charge species, and σ already includes their spin average and the identical scalar final factor. It was numerically checked in independent integration variables √s (outer) and dimensionless u=√s/T, plus angular 48/64/96/128 mapped refinements. For comparison assume illustrative RD H(T)=1.66 sqrt(g*ρ) T²/M_Pl (nonreduced M_Pl=1.2209e19 GeV) with **supplied** g* values, not derived SM Γ(T).

| T, illustrative g*ρ | Per-particle forward Γ_pair (GeV) | Γ_pair/H, mapped n128 |
|---|---:|---:|
| .1GeV, g*=17.25 | **3.916811688e−17** | **6.936017e3** |
| 1GeV, g*=61.75 | **1.030803042e−15** | **9.647821e2** |
| 100GeV, g*=106.75 | **3.287529814e−13** | **2.340226e1** |

For T=1GeV, an initial, **underresolved cosθ** numerical estimate was ~952H; properly mapped angular convergence gives **964.782H**. At T=100GeV, the old raw-angle calculation underestimated forward Γ/H by about a factor two, now **23.4023H**. At all stated Ts with y=.001, **the physically allowed pair channel's forward rate is not parametrically negligible relative to the illustrative Hubble scale** in the chosen MB test bath.

The exact vacuum cross section scales y⁴; at T1GeV, taking y=1e−9 instead of 1e−3 multiplies Γ_pair/H by 1e−24, leaving ~9.65e−22, so the old **small-κ weak-coupling FRW pilot** and the new strong-y benchmark can belong to qualitatively different 2↔2 regimes. The [independent Yukawa/Hubble comparison](PHI-BOSE-YUKAWA-HUBBLE-COUPLING-CONSISTENCY-2026-11-24.md) already showed the y=.001 Γrest/H(T=m1)≈1e11 versus the old adjustable κ≤.2. These are mutually consistent but do **not** imply physical DM abundance.

**Critical scope:** Γ_pair/H above is a **forward-only per-χ2 rate in an ideal MB bath with φ empty**. It is NOT the NET pair chemical-equilibration rate once inverse φφ→χ2antiχ2, Fermi statistics, (1+bφ)² and Bose occupations, evolving chemical potentials, SM thermal masses and portal interactions are included. Even Γforward/H≫1 by itself does not prove true chemical equilibrium, nor can a vacuum tree rate be naively extrapolated to the 100-GeV plasma without self energies/in-medium thresholds. Nevertheless, one cannot declare order-y⁴ pair processes unimportant for a full thermal model on grounds of perturbative order alone, when their forward rate in this specified benchmark is O(10³) H.

## Conservation and physical model consequences

The spin/channel construction proves that one can have χ1+antiχ1↔φφ and χ2+antiχ2↔φφ with the same allowed Yukawa interaction. The [integer 4×5 reaction rank certificate](PHI-BOSE-YUKAWA-CONTINUOUS-CHARGE-VS-ACCIDENTAL-COLLISION-2026-11-24.md) shows that these channels destroy the two *positive-count* collision invariants of the isolated 1↔2 solver; **only net Dirac U(1)** remains as a continuous linear charge in that set, while a possible dark Z2 parity is conserved modulo two. The existing 1↔2-only numerical conservation proofs remain valid *for that restricted model*; they are NOT a proof of two fundamental Noether charges in the full thermal plasma.

This creates a concrete model-design fork:
- If a genuinely tiny y is chosen to match the old weak-κ FRW benchmark, pair reactions may be negligible, but vacuum 1↔2 and nonlinear amplification must be re-solved at that y with actual SM portal.
- If the current y=.001 strong benchmark is retained, at least the allowed pair 2↔2 network must be **quantitatively** incorporated or excluded under the relevant in-medium thermal conditions, with correct energy and entropy. Actual number-changing channels can alter absolute Ω and the cold/hot fraction.
  
## Reproduction and remaining scientific gates

- [Independent covariant trace, mapped-angle σ and MB rate source](phi_bose_pair_annihilation_mb.py).
- [Direct 4×4 Dirac gamma-matrix trace, angular convergence and independently transformed MB rate tests](../tests/test_phi_bose_pair_annihilation_mb.py).
- Commands: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_pair_annihilation_mb.py; python -m kinetics.phi_bose_pair_annihilation_mb.
- Local independent scipy integration of mapped angle u used 48/64/96/128 orders: Γ_pair/H(T=1GeV) 965.1721/964.7965/964.78207/964.7820535; T100GeV 23.4052/23.40237/23.4022582/23.4022580. The plain-angle negative at √s100GeV is independently verified. Public Scientific CI must be checked for the corrected code/test commit.

This result still does NOT implement full quantum 2↔2 Boltzmann collision integral, its angle/momentum grid, SM reheating/thermal widths, a self-consistent expanding-background conserved spectrum, a 86-eV relic Ω, a physical condensate or CLASS likelihood.

**NEXT_ACTION:** derive the full microreversible φφ↔χ2antiχ2 and χ1antiχ1 pair collision kernels (including Pauli and Bose factors), thermal mass/in-medium propagators and any allowed χ–SM number-changing portals, then embed them alongside the existing physically normalized 1↔2 implicit scheme with **only truly conserved Noether charges**, actual energy/H and entropy. Quantify IR/UV/time collision convergence and compute real FRW fφ(q) before reusing the pre-existing verified CLASS thermal observations benchmark.
