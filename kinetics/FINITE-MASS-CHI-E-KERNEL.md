# Massive-electron invariant chi-e kernel: conditional tree-level extension (2026-11-21)

## Pre-registered scope and checks
Goal: replace the unphysical massless-electron threshold in the previously validated energy-transfer kernel without changing the vector-contact model. Independent limits: the old massless implementation; reversal of each physical four-body event; common-T Pauli detailed balance; forbidden transfer. Required numerical tolerances: massless rate relative 2e-13, reverse relative 2e-12, common-T Pauli relative 1e-13. **This is not an independent full 5D massive-electron normalisation or plasma heat validation.**

With M=mχ, m=me, p=√(E²−M²), k=√(ω²−m²), p'=√(E'²−M²), k'=√(ω'²−m²), and E+ω=E'+ω', the transfer magnitude ℓ obeys

    max(|p−p'|,|k−k'|) < ℓ < min(p+p',k+k').

Cone cosines for q=p−p'=k'−k:

    cχ=(p²+ℓ²−p'²)/(2pℓ),
    ce=(k'²−k²−ℓ²)/(2kℓ),
    t=(E−E')²−ℓ²,
    D0=2(Eω−pk cχ ce),
    <D²>=D0²+2p²k²(1−cχ²)(1−ce²).

For a parity-conserving four-Dirac-vector contact G(χbar γμ χ)(ebar γμ e), the spin-averaged tree-level trace gives

    |M|²=2G²[(s−M²−m²)²+(u−M²−m²)²+2t(M²+m²)]
         =2G²[2D²+2Dt+t²+2t(M²+m²)], D=s−M²−m².

After angular convolution, including massive energy-shell Jacobians d³p/(2E)→p dE/2, the same exact energy differential prefactor survives:

    K(E,ω;E',ω') = gχ ge/(256π⁵) ∫_(ℓmin)^(ℓmax) dℓ <|M|²>azimuth.

The contact model assumes a pointlike mediator and a **broken-electroweak Dirac electron current**. It is not an SU(2)L×U(1)Y-invariant high-temperature model by itself. Below/near MeV, electron/positron chemistry, photon reheating, annihilation, charge neutrality and temperature-dependent g* must be included before any thermal rate or FRW claim. An electron thermal population does not remain massless at T≪me.

## Reproduction and limits
Source: [finite_mass_chie_kernel.py](finite_mass_chie_kernel.py).
Regression: [test_finite_mass_chie_kernel.py](../tests/test_finite_mass_chie_kernel.py).
Run: `python -m pytest -q tests/test_finite_mass_chie_kernel.py` from repository root.

A standalone local implementation passed the massless-continuity, reversal, common-temperature Pauli, and forbidden-transfer controls before commit; the repository regression against the **prior independent massless code** awaits the new real GitHub CI. No measured massive-electron 5D heat benchmark or cutoff/grid convergence exists yet; no quantitative finite-me Q or γ_el is claimed. Prior energy-collocation h nonmonotonicity remains a negative result.

Next: check CI and correct any regression; independently integrate massive-electron CM 5D scattering heat at T/me={0.1,1,10}, compare the exact ℓ kernel with 3+ energy quadratures, and then derive electron/positron chemical and finite SM heat ledgers. In the early-universe EW-symmetric regime replace this broken-phase contact with an explicit gauge-consistent mediator or EFT matching. No 86-eV relic, BEC, ΩDM or CLASS-ready claim.
