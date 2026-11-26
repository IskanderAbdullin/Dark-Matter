# Globally correct **two-sided** massive Yukawa daughter-parent support, including 86-eV rest threshold (2026-11-24 UTC)

## Physical question

Earlier Yukawa boosted-parent weak tests use P=2/8/16 GeV and infer a scalar daughter minimum energy from the anti-parallel **endpoint** p=P. Within that restricted parent range the scalar's minimum energy decreases monotonically with p. This fails *globally* if the parent cutoff exceeds the boost needed to produce the daughter **at rest**: the anti-parallel minimum begins increasing at larger p, and a fixed energy is then supported only by a **finite** interval of parent momenta, with BOTH lower and upper physical bounds.

The previous conditional scalar-mode function was **only valid for P below the scalar rest-production parent threshold**. We derive and implement the globally correct interval, including χ2 and φ as two specializations. This is strict Lorentz kinematics and an independent negative check, not evidence of physical Bose condensation or relic density.

## General invariant-rapidity theorem

For the on-shell parent χ1 of mass M=m1=.1 GeV decaying isotropically to daughter species d of mass m_d>0 and rest-frame daughter momentum q*=KSTAR≈.048 GeV, define

    ξ_d = asinh(q*/m_d),
    ζ_d = asinh(sqrt(E_d²−m_d²)/m_d),
    y = asinh(p_parent/M), y_max=asinh(P/M).

The exact boosted daughter energy range at parent rapidity y is

    m_d cosh|ξ_d−y| <= E_d <= m_d cosh(ξ_d+y).

Equivalently at fixed daughter energy E_d, the parent is physically allowed **iff**

    |ξ_d−ζ_d| <= y <= min(ξ_d+ζ_d,y_max).

The exact two-sided parent momentum interval is therefore

    p_low=M sinh|ξ_d−ζ_d|,
    p_high=M sinh[min(ξ_d+ζ_d,y_max)],

provided p_low<=P. This is a *finite upper bound* whenever ξ+ζ < y_max. No event should be deposited from p outside that interval even if occupations are large; a sparse sampled graph may accidentally violate it. The formula is valid for both massive daughter χ2 and extremely light but positive-mass φ, independent of thermal FD/BE factors.

The **global minimum** daughter energy over parent p∈[0,P] follows from rapidity geometry:

    E_min,d(P)=m_d cosh(max(0,ξ_d−asinh(P/M))).

The daughter is producible **at rest** at the single parent momentum p*=M q*/m_d, and E_min=m_d for every P>=p*. If P>p*, using the anti-parallel energy evaluated at p=P as the **global floor is wrong**.

## Physical results, including a new prior-code counterexample

- For χ2 mass m2=.02 GeV, p*(χ2 at rest)=**0.23999999999979973 GeV**. At χ2 E=.0201 GeV with P2, the allowed parent range is **[.2151675203,.2672324797] GeV**; at χ2 rest-frame E*=.052 GeV it is **[~0,1.248] GeV**. A one-sided p_th..P generator would overcount unphysical states.
- For mφ=86 eV=86e−9 GeV, the parent boost needed to produce a scalar **exactly at rest** is **55,813.9534883 GeV**, about **55.8 TeV**. Thus all tested parent cutoffs P=2,8,16 GeV categorically exclude the physical φ zero-momentum rest region in this isolated decay channel.
- Physical scalar global energy floor for P=16 is **0.000149998547511642 GeV** (≈150 keV), not 86 eV.
- If P=60,000 GeV, scalar Eφ=86.1 eV is kinematically supported only by parent **p∈[53,186.4833,58,571.2236] GeV**, which ends *below* the 60,000 GeV cutoff. The old one-sided p∈[p_min,60,000] prescription would include forbidden ultra-boosted parent momenta.
- If P=100,000 GeV, the scalar anti-parallel energy at the upper parent endpoint rises to **1.01041666667e−7 GeV**, strictly ABOVE its true globally reachable rest mass 8.6e−8 GeV. The old practice kmin(P)=global minimum would falsely exclude a genuine on-shell daughter region. This does not mean that thermally populated GeV-scale parents are plentiful at ~55 TeV: their thermal occupation is exponentially tiny at T~1 GeV and no physical condensate inference follows from geometric accessibility.

## Implementation and independent checks

The new [globally exact kinematics module](phi_bose_invariant_parent_support.py) provides parent_interval(E,m,P) and global_support_floor(m,P). Its [three tests](../tests/test_phi_bose_invariant_parent_support.py) compare the rapidity interval with **independently reconstructed Lorentz-boosted daughter energy inequalities** for both daughters at P2/8/16; check χ2 finite upper bounds and the 55.8 TeV scalar threshold and two-sided range at P60TeV.

**Concrete corrective change:** updated [scalar fixed-bath modal generator](phi_bose_scalar_fixed_bath_strong.py) to integrate p across the complete two-sided parent interval, not p_low..P, and to use the globally correct scalar energy support floor. Original P2/8/16 results are unaffected. An initial overly strict equality test at exactly kmin(P) failed because rapidity inversion yielded y_low just above y_max by floating roundoff. A narrow 2e−13 *rapidity* tolerance now treats only such exact-boundary equality as a zero-width interval, without adding physically forbidden finite intervals. Eight combined local tests (the prior conditional scalar source tests plus these general kinematic checks) then passed in 3.05 seconds. Public CI separately verify.

Reproduce with OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_invariant_parent_support.py tests/test_phi_bose_scalar_fixed_bath_strong.py and python -m kinetics.phi_bose_invariant_parent_support. For P beyond 55 TeV this stage **verifies geometry only**: no ultra-high-energy finite thermal rate convergence, ultraviolet completion or cosmological extrapolation is claimed. At such cutoffs the previous parent-first cut-cell weak quadrature also needs a globally correct support partition, not just this scalar-mode correction.

## NEXT_ACTION / limitations

Continue the **fully coupled physically normalized** χ1/χ2/φ on-shell collision problem with actual three-species mode capacities, a common reversible event kernel, charge+energy and quantum entropy proof and finite-time Pauli/Bose positivity. For any UV parent extension beyond scalar p*, apply the two-sided rapidity theorem globally before integrating weak or strong sources. Prove or bound actual high-parent thermal tails; do not mistake formal IR accessibility for physically significant Bose amplification. No finite SM/FRW, dark relic Ω, CLASS or genuine condensate established here.
