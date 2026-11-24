# Gauge-consistent finite-vector electron portal and mandatory visible channels (2026-11-24 UTC)

## Scope and rationale

The connected Yukawa χ1 ↔ χ2 + real φ (mχ1=.1GeV, mχ2=.02GeV, mφ=86eV) admits a vectorlike dark U1 with q1=q2, NOT a diagonal χ1-only vector mediator with a neutral real φ. We examine one **conditional** broken-electroweak kinetically mixed massive vector A', with gD=.1, epsilon=1e−5 and mV=10GeV. It couples to **both** Dirac particle/antiparticle pairs and to the electromagnetic charged-lepton current. This is not a fully gauge/experimentally constrained UV model; the dark vector mass Higgs sector, electroweak/thermal mixing, running alpha, hadron spectrum and nonthermal SM are not derived. No physical 86-eV cosmic relic or CLASS-ready state follows.

## Distinct scattering and chemical channels

The actual spacelike t-channel unpolarized χ_i e− scattering amplitude uses the full physical factor (gD epsilon e)^2/(mV²−t)^2 and the exact Dirac trace. A 5D physical FD/Pauli weak heat kernel includes incoming lab momenta, relative angle, boosted CM polar and azimuth, fully blocked final χ/e and time-reversed events. All four separately counted Dirac charge sectors χ1,antiχ1,χ2,antiχ2 exchange energy with e±, with exact source Qdark+Qe=0, no elastic number creation and nonnegative microevent entropy. This is a proper instantaneous physical heat weak source for prescribed thermal FD distributions, **not** a nonlinear strong SM+dark finite-time relaxation operator.

At mχ2=.02 GeV,Tχ=.8GeV,Te=1GeV,mV=10GeV, the n=8,10,12,16 node-per-dimension five-dimensional results for one χ2 charge sector are Qχ2=(7.25883159,7.16316547,7.17093305,7.17087247)×10^−20GeV5. Last two relative difference 8.45×10^-6; n8 was +1.23% higher than n16. An *independent* mV→infinity contact source at n16 yields 1.135444594e−19GeV5, so finite mV suppresses heat to **0.6315475** of contact: a 36.85% difference. The physical full-vector spin trace also agrees with an independent complex 4×4 Dirac gamma trace. At a different four-sector benchmark Tχ1=.8,Tχ2=.9,Te=1,n8, total Qdark=+2.5872112410e−19 and Qe=−Qdark, H=+4.9034336263e−20 (appropriate GeV units). Older independent 5D finite-mediator studies already exist; this stage explicitly connects the vector portal to **both mandatory fermions**.

A DIFFERENT s-channel annihilation is χ_i+antiχ_i ↔ e−+e+ via the SAME massive vector. Its full forward MB rate is calculated with an actual finite ΓV including both χ1 and χ2 invisible partial widths, and an exact resonance transformation s=mV²+mVΓV tanθ. At T=1 GeV, χ2:
- Minimal (two-dark-plus-electron) ΓV = **5.3051646127×10^-3GeV**.
- Forward χ2 pair→e pair MB reaction density γ=**3.450407733×10^-16GeV4**; the corresponding forward-per-χ/H in a prescribed RD T,g* bath is 1593.812. This is **NOT NET chemical equilibration**; at equal temperatures and compatible chemical potentials reverse exactly cancels forward.
- Naive contact s-channel gives γ≈2.3972637314e−19, smaller by a factor **1439.31**, because its propagator misses the actual narrow on-shell vector pole.
- Independent spin-three real-vector narrow-width factorization n_A' Γ(A'→χχbar) Br(A'→e+e−) K1/K2 gives ≈3.447115846e−16 at T1, within 0.1% of the full Breit–Wigner integral; at T=.1GeV the same on-shell approximation fails because offshell production dominates. If A' is later a dynamical species, on-shell RIS cannot also be counted in the full 2↔2 kernel.

## Mandatory muon, tau and hadrons (negative completeness gate)

An electromagnetically mixed A' cannot consistently decay *only* to electrons when its mass is 10 GeV. At tree level for any charged lepton l:

    Gamma(A'→l+ l−) = alpha epsilon² mV/3
           * (1+2ml²/mV²) sqrt(1−4ml²/mV²),

with strict zero below threshold. At 10GeV e and mu are effectively massless; tau is open, so the full charged-lepton width is roughly **2.99× the electron-only width**. The corresponding on-shell NWA χχbar→all charged leptons also rises by about that factor (the total A' width remains almost entirely invisible for these declared gD/epsilon). Missing charged hadrons are more serious: the inclusive width must use a **measured R_had(mV²)=σ(e+e−→hadrons)/σ(e+e−→mu+mu−)** and the proper visible reference normalization, not free b-quark pole masses near physical bottomonium thresholds. The code therefore defaults to R_had=None and explicitly refuses to label its visible width complete. With an arbitrary illustrative parameter **R_had=3.5**, all charged visible widths would be about **6.5× the electron width**, but **3.5 has NOT been supplied as an actual measured value** and must not be used as a validated benchmark. Dark A' width changes negligibly because χ1 and χ2 dominate at these couplings. See dark-photon charged-lepton width and data-driven hadronic R treatment in [PTEP 2026 013C02](https://academic.oup.com/ptep/article/2026/1/013C02/8404152) and [Phys. Rev. C 112 054905 (2025)](https://journals.aps.org/prc/abstract/10.1103/qzzg-h1xw).

## Permanent code/evidence and further gates

- [Physical finite t-channel source](phi_bose_vector_portal_finite_mediator.py), [source/5D Pauli heat](phi_bose_vector_portal_pauli_heat.py), [resonance/inverse/reverse energy code](phi_bose_vector_portal_annihilation_resonance.py), [visible lepton/Rhad code](phi_bose_vector_portal_visible_width_gate.py).
- [Dirac amplitude/transport tests](../tests/test_phi_bose_vector_portal_finite_mediator.py), [Pauli 5D heat tests](../tests/test_phi_bose_vector_portal_pauli_heat.py), [resonant/NWA tests](../tests/test_phi_bose_vector_portal_annihilation_resonance.py), [visible threshold/Rhad tests](../tests/test_phi_bose_vector_portal_visible_width_gate.py).
- Independent Sci CI: five-field focused runs #45430331573 (Pauli), #45430926445 (refined heat), #45431671479 (resonant inverse/entropy), #45432164795 (NWA) and #45452679130 (visible) passed; whole dual-version Scientific status to be checked by SHA.
- [Machine-readable fixed benchmarks](results/phi-bose-vector-portal-partial-sm-2026-11-24.json). **No** physical finite electron mass at T~me, nonzero muon/hadron plasma heat, dynamical A' abundance, ΓV(T) thermal loop, BBN, direct detector/collider constraints, complete five-sector dark+SM time evolution, actual freezeout fφ(q), Omega or CLASS likelihood has passed.

**NEXT ACTION:** specify a true phenomenologically admissible vector/gauge completion including visible hadronic spectral input, dark Higgs/thermal width and relevant SM channels, then build a single reversible FD/BE+SM collision operator with E/Q/H and late FRW before projecting any chemical or kinetic Γ(T).
