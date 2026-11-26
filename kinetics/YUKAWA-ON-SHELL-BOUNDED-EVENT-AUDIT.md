# Bounded on-shell χ1↔χ2+real φ quantum event: exact energy, Pauli/Bose and entropy (2026-11-22)

## Scientific question and scope

The prior [moving on-shell Yukawa source](MOVING-YUKAWA-BOSE-PAULI-WEAK-AUDIT.md) validates boosted daughters and quantum number/heat **WEAK** rates, but no physical evolving φ distribution exists in the two-χ shared-lepton solver. Before attempting a full physical q grid, test a single on-shell microscopic χ1↔χ2+φ reaction as a strictly bounded quantum elementary event. The capacities of its three modes are **arbitrary positive test weights**, not a physical density-of-states mapping. This is a useful but incomplete strong-collision building block, not a cosmological rate.

Let E1=E2+Eφ, f1=a−x/w1, f2=d+x/w2, bφ=b0+x/wφ. Then the **exact microphysical** Pauli/Bose forward-inverse factor is

    J(x)=f1(1−f2)(1+bφ)−(1−f1)f2 bφ
        =f1(1−f2+bφ)−f2 bφ,

a quadratic polynomial in the weighted event extent x. For arbitrary **positive test reaction weight** W and timestep dt, define a **bounded first-order implicit** event by

    x=dt W J(x),
    x_lower=max[(a−1)w1, −d w2, −b0 wφ],
    x_upper=min[a w1,(1−d)w2].

Solve by bracketed Brent on the physical extent interval, without post-hoc clipping. This is **not** an exact continuous-time flow, global multi-event integrator, nor a proof of unique implicit root for all possible physical/stiff configurations.

Three exact linear invariants of this ISOLATED 1↔2 event are

    C_f=w1 f1+w2 f2,
    C_1φ=w1 f1+wφ bφ,
    E=w1 E1 f1+w2 E2 f2+wφ Eφ bφ.

The fermionic/Bose entropy is

    S_F=−w[f ln f+(1−f)ln(1−f)],
    S_B=wφ[(1+b)ln(1+b)−b ln b],

and its microscopic derivative follows W(F−R)ln(F/R)≥0. The finite entropy is well defined via xlogy even when the initial Bose occupancy is exactly zero; no infrared Bose condensate is inferred.

## Reproducible on-shell numerical experiment

Get three exactly energy-matched levels from the earlier **moving scalar Yukawa Lorentz boost** with m1=.1,m2=.02,mφ=.001,mother momentum p=.3 and daughter cosθ*=.2:

    E1=.31622776601683794,
    E2=.1932193769133251,
    Eφ=.12300838910351285.

Choose synthetic capacities (w1,w2,wφ)=(.1,.15,.2), positive coefficient W=.05 and dt=10. These capacities are not claimed to match an actual finite phase-space cell, so the numbers below are dimensional/methodology checks.

| Initial (f1,f2,bφ) | Implicit x | Final (f1,f2,bφ) | ΔS |
|---|---:|---|---:|
| (.6,.3,.1) | +.0328142120 | (.27185788,.51876141,.26407106) | +.0659934336 |
| (.05,.8,2) | **−.0365123342** | (.41512334,.55658444,1.81743833) | +.0605746330 |
| (.4,.1,0) | +.0284780807 | (.11521919,.28985387,.14239040) | +.0958931627 |

The inverse reaction really **destroys φ** when sufficiently occupied; the empty-boson reaction correctly creates φ. The two charges and total collocated energy differ by ~10⁻¹⁷. At common T=.8 and μ2=.02, μφ=−.01, μ1=.01=μ2+μφ, the on-shell thermal Fermi/Bose state is stationary to <1e−13.

The deterministic seed=20261122 stress test covers **250** internal Fermi pairs in [.001,.999], Bose b=10⁻⁵–10⁴, each positive capacity 10⁻³–10², reaction W 10⁻³–10² and dt 10⁻⁴–10³. No tested case violated Pauli/Bose, entropy, or the two charges and energy within 5e−14 relative plus 5e−13 absolute numerical precision. An initial overstrict **absolute-only** tolerance 3e−12 flagged a ~3.64e−12 rounding discrepancy on a quantity ~1.9×10⁴ (relative ~2e−16). Corrected the scientific test to use both relative and absolute tolerances; no physics criterion was relaxed. Five local tests passed.

## Reproduction, negative limitations and next gate

[Single-event source](yukawa_on_shell_bounded_event.py), [five tests](../tests/test_yukawa_on_shell_bounded_event.py), [numeric evidence](results/yukawa-on-shell-bounded-event-2026-11-22.json), [dated stage](../archive/stages/2026/11/2026-11-22-yukawa-bounded-on-shell-event.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_yukawa_on_shell_bounded_event.py

**Critical limitation:** proof of on-shell energy/entropy of one synthetic mode triplet **does not provide** physical daughter density-of-states capacities, angular projection, rate matching with continuous moving 1↔2, bounded global multi-event time stepping, physical φ zero-mode or actual abundance in the shared-lepton coupled plasma. A naïve energy barycentric remap can violate entropy, as demonstrated earlier in the χℓ sector. Next P0 is an energy-matched finite-volume χ1/χ2/φ strong event graph with physical weights and 3+ q/angle/dt/IR/cutoff checks, then couple with the two-χ shared lepton bath and SM/Zprime/FRW. No physical Γelastic, 86eV relic f(q), ΩDM, BEC or CLASS-ready result. Status in_progress.
