# Is the benchmark Yukawa y=10^-3 compatible with the weak-κ FRW pilot? (2026-11-24)

## Scientific question

The new three-species physical-capacity χ1↔χ2+φ continuum/finite-time code uses **actual on-shell tree Yukawa coupling y=0.001**, masses m1=.1GeV, m2=.02GeV, mφ=86e-9GeV. The older [FRW fixed-bath pilot](FRW-FIXED-BATH-PILOT.md), however, introduced an *independent adjustable* κ=Γrest/H(T=m1) and ran κ=1e-4,1e-3,.01,.1,.2. Does the concrete y=.001 benchmark actually inhabit that old weak-κ parameter range in a radiation-dominated universe?

**No, very decisively**, for the specified vacuum Lagrangian and standard illustrative g* background. The older pilot cannot be used as the cosmological evolution of the new y=.001 numerical collision model. This is a simple but important microscopic/FRW **consistency gate**, not a solved SM thermal history or exclusion of the whole model.

## Microscopic tree-level width

For a Dirac χ1→Dirac χ2+real scalar φ with interaction L_int=−yφ χbar2 χ1+h.c., unpolarized spin average yields

    |M|²_averaged = y²[(m1+m2)²−mφ²],
    kstar = sqrt(λ(m1²,m2²,mφ²))/(2m1),
    Γrest = |M|²_averaged*kstar/(8πm1²)
          = y²[(m1+m2)²−mφ²]sqrt(λ)/(16π m1³).

This matches the width used in the independently verified Yukawa weak physical source. At y=.001:

- Γrest=**2.75019741662424e−9 GeV**;
- vacuum rest lifetime ℏ/Γ=**2.3933262133e−16 seconds**.

The width is a **vacuum tree-rate**, not an actual thermal plasma Γ(T) or an SM-equilibration rate; the illustrative g1=g2=4 counts a Dirac fermion and antifermion equally populated, with gφ=1 one real scalar. A full spin-channel/portal theory must separately justify that convention.

## Standard RD Hubble comparison and decisive parameter mismatch

Under a specified radiation-dominated background with **constant supplied** g*ρ(T),

    H(T) = 1.66 sqrt(g*ρ) T² / M_Pl,    M_Pl=1.2209e19 GeV
(non-reduced Planck mass). The *same* old pilot definition is κ(T)=Γrest/H(T). Using g*ρ=17.25 for an illustrative T=m1=.1 GeV gives H=5.64706184e−21 GeV and

    κ(T=m1, y=.001)=**4.8701386531e11**,

whereas the old FRW pilot ranged from κ=.0001 to .2. Even if one varies the uncertain g*ρ(.1GeV) from 17.25 to 61.75, κ remains **2.5740525e11 to 4.8701387e11**; this is many orders of magnitude beyond its tested region. To obtain the old κ=.1 at T=.1GeV and g*ρ=17.25 with the same vacuum masses, one would need

    y_required = .001 sqrt(0.1 / 4.870138653e11)
               = **4.53136804e−10**,

over 2.2 million times smaller in amplitude coupling. For g*=61.75 the same target y≈6.23292e−10. One CANNOT both retain y=.001 in the physical collision source and silently reuse a κ=.1 FRW pilot result for absolute Ω/freezing.

At T=1GeV,g*=61.75, H=1.06843099e−18GeV, Γrest/H≈2.574052465e9 (also far outside κ≤.2). At g*=106.75 higher temperatures the ratios decrease strongly because H∝T².

## Critical relativistic time-dilation caveat

The bare ratio Γrest/H is **not** itself the thermal per-parent decay-rate/H for relativistic χ1. For dilute, ideal Maxwell-Boltzmann thermal parents, an analytically exact independent average is

    <m1/E1>_MB = K1(m1/T)/K2(m1/T),
    <Γ_lab>_MB = Γrest*K1/K2.

The Bessel identity was independently checked by positive thermal p²dp integrals of exp(−E/T) and (m1/E)exp(−E/T) for T=.1,1,100GeV (test tolerance 2e−10). For T≫m1, K1/K2≈m1/(2T). This test-parent approximation ignores quantum fermionic blocking, stimulated Bose inverse processes, actual χ2/φ distributions, SM heat, thermal masses and changing g*S(T).

| Temperature, illustrative g*ρ | Γrest/H | <m1/E1>_MB | <Γ_lab>_MB/H |
|---|---:|---:|---:|
| T=.1GeV, g*=17.25 | 4.87014e11 | .37044117 | **1.80410e11** |
| T=1GeV, g*=61.75 | 2.57405e9 | .04939172 | **1.27137e8** |
| T=100GeV, g*=106.75 | 1.95773e5 | .000499998 | **97.886** |
| T=500GeV, g*=106.75 | 7.83090e3 | .00009999998 | **0.78309** |
| T=1000GeV, g*=106.75 | 1.95773e3 | .000049999998 | **0.097886** |

For a **fixed** g*=106.75 RD bath, y=.001 and the MB *test-parent average only*, the equation <Γ_lab>_MB/H=1 occurs at **T≈460.865GeV**, NOT at the many-TeV Γrest/H=1 crossover. This is a diagnostic, **not** a proved real thermalization/freeze-out temperature; actual inverse reactions, density, chemical affinities and portal Γ(T) can profoundly change the collision dynamics. For T=.1 or 1GeV, even including time dilation leaves the illustrative per-parent decay timescale far faster than H, emphasizing that weak open-source cosmological freeze-in pilot parameters and this dark-sector collision benchmark are distinct hypotheses.

## Link to 86-eV IR support and cosmological caveat

The earlier rigorously derived p≈55.814TeV parent momentum for producing exactly-at-rest 86-eV scalar in THIS 1↔2 channel is a **kinematic** number, not a production epoch or a parent temperature. The rigorous exponentially tiny **fixed** T1=.9GeV FD/BE near-rest rate bound [86–86.1eV rate <10^-23093GeV4](PHI-BOSE-NEAR-REST-THERMAL-ABSOLUTE-BOUND-2026-11-24.md) therefore **does not apply to an early hypothetical T~tens of TeV plasma**: the parent tail, H, Lorentz time dilation and occupations are all different there. Nor does the MB crossover T~461GeV alone imply a physical zero-mode occupation or thermalized 86-eV relic. A genuine early-Universe result requires full χ+SM portal rates, number changing and finite SM entropy/energy with FRW.

The two collision invariants N1+N2 and N1+Nφ of the restricted **single 1↔2 event network** are exact for that reaction, but they are NOT automatically two independent fundamental global Noether charges of the full Yukawa field theory: adding physical 2↔2/annihilation/SM channels may break the second invariant. In full cosmological model explicitly rederive which charges are truly conserved and which are only accidental stoichiometric invariants.

## Reproduction and status

- Source: [phi_bose_yukawa_hubble_gate.py](phi_bose_yukawa_hubble_gate.py); four focused independent tests [test_phi_bose_yukawa_hubble_gate.py](../tests/test_phi_bose_yukawa_hubble_gate.py) reproduce Dirac spin trace vs old width, Maxwell-Boltzmann K1/K2 vs independent direct physical momentum integrals, κ mismatch for g* band, and high-T crossover.
- Run commands: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_yukawa_hubble_gate.py; OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_yukawa_hubble_gate.
- Scientific CI of new source/test commits must be checked separately. This is a *comparison to supplied RD parameters*, not independent validation of experimental laboratory Yukawa model or an in-medium width.

**NEXT_ACTION:** choose/derive the actual microscopic benchmark parameters and physical χ–SM portal with non-arbitrary Γ(T), real g*(T), fast/slow channel assessment and full energy-number exchange; then use the *existing* coupled physical-capacity implicit dark-sector solver as an interaction substep of a genuinely closed FRW evolution. Keep class_ready=false and do not derive real Ω/BEC/CLASS from a free old κ pilot.
