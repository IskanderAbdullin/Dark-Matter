# EXACT real-scalar massive 2↔2 Yukawa pair IR support: χ2-pair reaches 86eV rest at 4.65TeV parents, not 55.8TeV (2026-11-24)

## New on-shell kinematic result and relation to older 1↔2

The original χ1→χ2+realφ vacuum two-body process with m1=.1GeV,m2=.02GeV,mφ=86e−9GeV reaches φ rest only for parent momentum p≈55.81395TeV. We have now derived a separate **globally exact support theorem** for allowed **χi+antiχi→φ+φ** on-shell pair collisions, when each incoming equal-mass Dirac particle has hard momentum cutoff p≤P and mχ>mφ. The IR support can change dramatically even at modest P.

For two incoming equal-mass fermions of mχ and bounded individual p≤P, the total initial E obeys E≤2sqrt(P²+mχ²), while CM invariant s≥4mχ². In this LAB frame, a scalar with mass mφ among the two final identical scalars has minimum energy at fixed E and s:

    k_min(E,s) = [E - sqrt(E²-s)*sqrt(1-4mφ²/s)]/2.

The boost term is maximized at invariant s=2mφ E if that invariant is allowed, otherwise at the minimum s=4mχ². For P below the rest-access threshold, 2mφ E≤4mχ² always, so the exact global minimum occurs when both incoming parents have the same maximal momentum P and are perfectly collinear, hence CM relative momentum zero, with one outgoing scalar emitted backwards. This provides an ATTAINABLE exact minimum

    k_min_pair(P) = sqrt(P²+mχ²) - P sqrt[1-(mφ/mχ)²].

At large P direct subtraction loses precision. The numerically stable, **globally valid** rapidity expression is

    k_min_pair(P) = mφ cosh[max(0,
                             arcosh(mχ/mφ) - arsinh(P/mχ))].

Define p_star=mχ sqrt(mχ²−mφ²)/mφ ≈mχ²/mφ. At P≥p_star, one may choose both incoming parents at p=p_star and emit a scalar EXACTLY at rest, so the global lower limit clamps at mφ. The endpoint evaluation with both parents p=P would misleadingly rise again beyond p_star; **not** the correct global support.

## Quantitative physical massive-eV results

For mφ=86eV, mχ2=.02GeV, mχ1=.1GeV:

| Each incoming parent cutoff P | χ2 antiχ2 pair minimum φ energy | χ1 antiχ1 pair minimum φ energy | Restricted original χ1 decay scalar floor |
|---:|---:|---:|---:|
| 2GeV | **99.9975186keV** | **2.49843945MeV** | ~1.199251MeV |
| 8GeV | **25.0000349keV** | **624.975591keV** | ~299.988keV |
| 16GeV | **12.5001430keV** | **312.496954keV** | ~149.998keV |

The *lighter* χ2-pair reaches much lower scalar energies than the heavier χ1 pair and the original single decay. At P2–16 its IR floor is ~12 times lower than decay. Exact scalar rest for the **χ2 pair** becomes kinematically possible if each incoming χ2 has p≥**4651.16279065GeV** (~4.65TeV), whereas χ1 pair needs p≥**116279.069767GeV** (~116TeV per parent); original 1↔2 needs ~55.8TeV for the parent. This is a **kinematical accessibility fact**, not a reaction rate, thermally populated 86eV condensate or proof of cosmological cold DM.

At T=1GeV fixed thermal FD distributions, pair-rest source parents p≥4.65TeV per particle are exponentially suppressed by at least exp[-2×4651] in their combined MB Boltzmann factor, so an accessible kinematic shell is NOT necessarily efficiently populated. At hypothetical early T~TeV this simple suppression can be much less severe, but actual rates then require consistent in-medium masses, elastic widths, SM entropy and FRW; no extrapolation of fixed-GeV IR rate bounds is justified.

## Crucial finite-volume consequence for the existing five-field strong solver

The first five-species coupled Bernstein-P2 prototype initially reused scalar grid cut-cells based only on **χ1 1↔2** parent support. This leaves a broad φ cell from mφ to the 1↔2 floor (≈1.199MeV at P2), which is **NOT entirely collision-inactive** once χ2 pair is included: actual 2↔2 physical daughter events begin near **100keV**. Without inserting the pair-specific support boundary the quadratic physical Bose entropy-dual projection can smear the active source into the truly unpopulated lowest interval and create spurious IR inverse-η tangents. The support-correct grid must include BOTH pair floors plus decay floor as exact φ energy edges. This is a structural correction, not a cosmetic resolution change. It also leaves modes below all physical supports explicitly collision-frozen under the hard P.

This is consistent with earlier proven non-positive ηdot polynomial projection counterexamples: the positivity of Bernstein hats and weak event conservation do not imply globally sign-preserving strong inverse-fugacity tangents or macroscopic BEC.

## Proof, precision, and limitations

Module [positive massive pair IR floor](phi_bose_pair_ir_parent_support.py) implements direct analytical and stable rapidity support formula and exact parent-rest momentum. Four focused tests [independently minimize over physically conservative bounds](../tests/test_phi_bose_pair_ir_parent_support.py) across P2/8/16 and two χ masses, verify global lower bound from (E,s) minimization, and provide **above-p_star false-endpoint negative control**. Primary numerical table [JSON](results/phi-bose-pair-ir-support-2026-11-24.json) and [dated archive](../archive/stages/2026/11/2026-11-24-pair-massive-scalar-IR-support.md). For very small mφ/P use rapidity, not double-precision cancellation.

Important boundary: proof assumes exactly equal-mass incoming χi and antiχi, isotropic on-shell **vacuum** 2→2, both individually truncated by the same hard momentum P; masses satisfy mχ>mφ. It states only support of the full phase space; an event quadrature with few angle nodes will not hit its exact mathematical lower endpoint. Finite-T/SM pole masses and physical rates modify the relevant support. A completed late-time quantum/FRW model with actual mass spectrum is not yet available.

NEXT_ACTION: update the five-species real-Bose physical-capacity grid to split at the pair χ2 and χ1 physical shell floors in addition to decay; verify frozen below-support scalar modes, finite-time positive Bose realizability, independent q/angular/time grid convergence and actual Γ/H. Later derive in-medium χ/φ, SM bath, Friedmann and absolute frozen f(q), Ω and observational likelihood. No CLASS-ready/BEC claim.
