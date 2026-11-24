# TWO dark numerical number charges are not TWO fundamental Yukawa charges (2026-11-24)

## Why this matters physically

The physically normalized χ1↔χ2+φ three-species continuum and implicit Bernstein-P2 solvers exactly conserve **Nχ1+Nχ2** and **Nχ1+Nφ** (with gχ1=gχ2=4 lumping particles and antiparticles). These are valid stoichiometric invariants of the *isolated unpolarized single 1↔2 channel*. They do NOT automatically imply two independent Noether charges of the actual microscopic dark matter model.

Exhibit a concrete, minimal illustrative renormalizable theory:

    L_int = −y φ bar(χ2) χ1 + h.c.

with Dirac χ1,χ2 (masses .1,.02GeV), one REAL φ (86eV), and y=.001. A possible dark Z2 assigns χ1 and φ odd, χ2 even, protecting φ stability as the lightest odd state if other interactions honor Z2 and no parity-breaking vacuum occurs. Both Dirac fermions share an independent vector-like U(1) charge +1, anti-fermions −1, real φ charge0. An actual portal, gauge anomaly and self-interaction specification are still absent; this is only an illustrative allowed minimal symmetry structure.

## Exact alternative reaction stoichiometry

Keep five distinguishable populations (χ1,antiχ1,χ2,antiχ2,φ). The same Yukawa vertex permits the two inverse-decay reactions plus 2↔2 scalar-pair reactions via virtual **t/u-channel fermion** exchange. Do NOT assume their amplitudes are thermal-rate dominant; they arise at order y^4 and need a separate actual squared-amplitude and finite-temperature comparison, including thresholds and Pauli/Bose factors.

| Channel | Net stoichiometry in ordered five species |
|---|---|
| χ1→χ2+φ | (−1,0,+1,0,+1) |
| antiχ1→antiχ2+φ | (0,−1,0,+1,+1) |
| χ1+antiχ1→φ+φ (virtual χ2) | (−1,−1,0,0,+2) |
| χ2+antiχ2→φ+φ (virtual χ1) | (0,0,−1,−1,+2) |

The restricted two-row decay stoichiometric matrix has rank **2**, nullity **3** in these five species (one physical vector U(1) and two accidental additional linear counts for this restricted network). Adding the two pair channels yields a 4×5 integer S with exact minor det[S_columns(χ1,antiχ1,χ2,φ)]=**−2**, hence rank4 and nullity1.

The unique continuous linear conserved charge (up to scale) for these four channels is

    Q_U1 = Nχ1 − Nantiχ1 + Nχ2 − Nantiχ2,

because S·(1,−1,1,−1,0)^T=0 exactly. Candidate discrete Z2 parity Nχ1+Nantiχ1+Nφ modulo 2 also survives, but it is **NOT** another continuous particle-number charge.

The numerical old positive total Dirac counts and old Nχ1+Nφ extended to particle+antiparticle counts are q_F=(1,1,1,1,0), q_1φ=(1,1,0,0,1). They satisfy decay-only S0·q=0, but the two pair channels give

    S_full q_F  = (0,0,−2,−2),
    S_full q_1φ = (0,0, 0,+2).

Both previously reported positive-count invariants can fail in a physically allowed wider network. This is a mathematically exact **scoped model-closure negative result**, not proof that pair channels are cosmologically important for y=.001 or in a more general portal theory.

## Consequences for dark-cosmology research

1. The current fully coupled 1↔2 dark-sector implicit solver remains **correct for the model that was actually discretized**: two stoichiometric invariants and energy at each on-shell event. Do not weaken those source/unit tests.
2. A realistic relativistic number/entropy/energy budget for φ dark matter may be modified by 2↔2 pair creation/annihilation. This can change whether the eventual absolute Ω can be inferred from initial fermion inventory. The published 2026-11-21 FRW fixed-bath pilot and old reservoir bound are valid only under their specified restricted 1↔2 assumptions.
3. Quantitatively derive actual vacuum/thermal amplitudes and 2↔2 collision rates, check required SM channels and Z2/U1 assignments, and show when/if these channels may be discarded. Symbolic y^4 power counting is NOT sufficient because H and particle densities depend on T, thresholds and effective thermal masses. The [Yin 2023 1↔2 cold production paper](https://doi.org/10.1007/JHEP05(2023)180) and [Sakurai–Yin 2025 stimulated 2↔2 scattering](https://doi.org/10.1007/JHEP03(2025)202) are physical motivations, NOT proofs of this particular pair-annihilation rate.
4. Once the SM portal and kinetic equations are defined, use the exact continuous symmetry (net Dirac number) and discrete Z2, not simply both old positive-number invariants, to validate a fully open physical cosmological collision network. Entropy and total energy still need the actual SM/FRW heat channels.

## Reproduction, scientific status

- [Exact 4×5 integer stoichiometric module](phi_bose_yukawa_charge_audit.py); [four independent invariant/nullspace tests](../tests/test_phi_bose_yukawa_charge_audit.py).
- Commands: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_yukawa_charge_audit.py; python -m kinetics.phi_bose_yukawa_charge_audit.
- The rank certificate includes the explicit nonzero integer determinant −2 and S·q_U1=0. It does not numerically compute an actual amplitude for a forbidden-by-chirality case, the 2↔2 rate, or a full U(1) gauge/portal completion. Verify public Scientific CI separately.

**NEXT_ACTION:** derive and test the physical χ1/χ2/φ dark+SM Lagrangian (continuous charges, Z2 stability, thermal masses/mediator) and at least one pair 2↔2 amplitude/thermal rate plus the original 1↔2. Quantify when 1↔2-only stoichiometric charges approximate the full kinetics before declaring a stable physical 86-eV relic and passing it to the already verified thermal CLASS pipeline.
