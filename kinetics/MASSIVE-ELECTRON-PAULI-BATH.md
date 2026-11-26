# Physical electron mass in the direct microreversible full-Pauli χe kernel

**Scope:** controlled Dirac χ mass M and physical electron/positron mass m_e=0.00051099895 GeV in the contact-vector χe 2→2 Lorentz-invariant collision. This is a conditional massless-to-massive e± submodel, **not** real early-Universe SM plasma, 86-eV Bose DM, BEC or CLASS.

## Independent trace and invariant angular derivation

The unequal-mass spin-averaged vector amplitude is

    |M|² = 2 G²[(s−M²−m_e²)²+(u−M²−m_e²)²+2t(M²+m_e²)],
    s+t+u=2(M²+m_e²).

For lepton energies ε,ε′ with momenta k=√(ε²−m_e²), k′=√(ε′²−m_e²), and χ momenta p,p′, let ℓ be the three-momentum transfer, t=(E−E′)²−ℓ², and
cχ=(p²+ℓ²−p′²)/(2pℓ), ce=(k′²−k²−ℓ²)/(2kℓ).
Define D0=2(Eε−pk cχce), and azimuth average D²=D0²+2p²k²(1−cχ²)(1−ce²). Then

    <|M|²> = 2G²[2<D²>+2D0t+t²+2(M²+m_e²)t].
    K(E,ε;E′,ε′)=gχge/(256π⁵) ∫_(ℓmin)^(ℓmax) dℓ <|M|²>,
    ℓmin=max(|p−p′|,|k−k′|), ℓmax=min(p+p′,k+k′).

This follows from the same exact invariant four-body spherical shell convolution as the previously independently verified massless-electron q kernel. Energy is exactly conserved, E+ε=E′+ε′. The kernel is **directly** time-reversal symmetric without post-hoc rate arithmetic averaging. Source: [microreversible_massive_lepton.py](microreversible_massive_lepton.py); [tests](../tests/test_microreversible_massive_lepton.py). As m_e→0 it returns the earlier independent massless kernel to numerical roundoff.

Physical radial Fermi capacities for each charged lepton energy interval are ge(k_hi³−k_lo³)/(6π²), rather than k³ for a massless particle. A common energy lattice for χ and e still conserves each quartet centre's total energy exactly and obeys Pauli F/R. This is an energy-collocation approximation, not exact within-cell full Boltzmann.

## Tested charge-symmetric pair-bath finite-mass effect

For a deliberately specified FD test with Mχ=10 MeV, G=0.001 GeV⁻², Tχ=0.8 Te, μχ=Mχ−2Te, μe=0, gχ=ge=4, nχ=ne=100, energy spacing Te/4, and 8 transfer Gauss nodes, finite-mass and massless electron χ-heating Qχ ratios are:

| Te (MeV) | m_e/Te | Qχ(m_e)/Qχ(m_e=0) |
|---:|---:|---:|
| 2.0 | 0.2555 | ~0.9942 |
| 1.0 | 0.5110 | ~0.9857 |
| 0.5 | 1.022 | ~0.9651 |
| 0.2 | 2.555 | ~0.8097 |
| 0.1 | 5.110 | ~0.35563 |
| 0.05 | 10.220 | ~0.01968 |

At Te=.1 MeV: Qχ(m=0)=3.7892469043×10⁻43 GeV⁵ vs Qχ(m_e)=1.34758×10⁻43 GeV⁵. The independent electron/positron FD **number-density** suppression ratio there is only ~0.06755; the heat suppression is ~0.35563 because the heat-weighted collision kernel disproportionately samples energetic electrons. At Te=.05 MeV the electron number suppression is ~0.0009855 vs heat suppression ~0.01968. **Negative result:** replacing the massive-electron heat operator by the massless-electron operator multiplied by just n_e(m)/n_e(0) is quantitatively invalid for this specified model.

The Te=.1 MeV h=Te/4 mesh ratio at n=60,80,100,140 is ~0.359341, 0.355770, 0.355633, 0.355630: n100→140 is stable to ~1e−5 relative, but n60 is ~1% off. This is an explicitly limited numerical convergence check, not a cosmological inference.

Crucially this models an **ideal charge-symmetric e+e− bath at zero chemical potential only**. After cosmic positron annihilation, ordinary electrons associated with baryons **remain**; their density, µe and heat exchange are not represented here. Electron annihilation/photon entropy release, neutrino decoupling, weak chemical processing, QED plasma thermal masses, UV gauge completion and g*ρ/g*S/H are omitted. This scan has Mχ/Te=5…200, not the original relativistic BEST Mχ/Te=.1, and cannot calibrate that BGK. The absolute toy heat rates are not dark-matter relic abundances.

**Reproduce:** python kinetics/verify_massive_lepton_threshold.py --output outputs/massive-e-pauli.json. [Compact results](results/massive-e-pair-threshold-2026-11-21.json). The immediate next physics gate is joint e−/e+/ν lepton chemistry plus finite SM heat, complete gauge/mediator model and independent 5D collisional checks; only thereafter consider coupling to pinned BEST χ1↔χ2+86eV φ FRW freezeout.
