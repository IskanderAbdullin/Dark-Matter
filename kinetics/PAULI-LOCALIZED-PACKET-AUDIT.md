# Localized Pauli χe collision packets: a slower **variational** IR mode, and a dangerous quadrature trap

**Hypothesis:** Global polynomial perturbations do not necessarily find the slowest response of a physical fermionic phase-space distribution, especially near p≈0. Test sharply localized, number-orthogonal logit perturbations using the same exact 5D full-Pauli collision Dirichlet benchmark. Do **not** claim an eigenmode or a finite/zero *true* kinetic spectral gap from one trial packet.

## Test function and conservation

Let FD equilibrium have T=1 GeV and Mχ=.1 GeV, electron-only vector-contact G=.001 GeV⁻², gχ=ge=4. Test `h(E)=exp[-(E/T−center)^2/(2 width^2)]−<exp(...)>_FD_susceptibility` where the average uses the exact χ number susceptibility weight `p² fFD(1−fFD)`. The subtraction removes χ's conserved number null mode. The Rayleigh quotient

```text
R[h] = D[h,h]/S[h,h],
S[h,h]=gχ/(2π²)∫dp p² f(1−f)h(E)^2,
D[h,h]=(gχ ge/2)∫d³p d³k/(2π)^6 vM ∫dΩ* dσ/dΩ*
                 F_eq [h(E')−h(E)]²,
```

is positive and provides a conditional **variational upper bound** on the infimum of the stated symmetric kinetic collision operator if its domain assumptions hold. It is not a directly computed eigenvalue or a proof that a nonzero lower bound on the spectrum exists.

## Negative quadrature result and correction

A narrow IR packet centered at E/T=.25 with width=.08 mostly lives around p/T≈.229. A naive Gauss quadrature of just n=20 nodes spanning p/T∈[0,32] **misses** this physically important support: it yields R≈**2.60445×10⁻8 GeV**. Changing that uniform n to 16,24,32,40,48 makes the apparent Rayleigh quotient oscillate substantially (for example n32 ≈6.002×10⁻8; n48 ≈5.558×10⁻8). **That 2.604×10⁻8 value is a numerical artefact, not a reliable slow kinetic rate.**

Implement [localized_pauli_packets.py](localized_pauli_packets.py) with a composite Gauss scheme splitting **both incoming χ and electron momentum integrals at 2T**, retaining exact CM scattering and final-state Pauli factors. For the **same** packet and cutoff32:

| IR-adapted angular n, per momentum segment 2n | Rayleigh, GeV |
|---:|---:|
| 12 | 4.710313×10⁻8 |
| 16 | 4.857177×10⁻8 |
| 20 | 4.856938×10⁻8 |
| 24 | 4.853594×10⁻8 |

The difference n20→24 is ≈0.069% while cutoff20→40 at n20 changes by under 0.004%. The separate packet center E/T=.15,width=.15 gives R≈4.57518×10⁻8 (n20) and 4.56958×10⁻8 (n24) and is slower still. These are **empirically convergent numerical Rayleigh quotients for exactly specified trial functions**, not proof that the full operator possesses an eigenvalue at those rates or a gap of that size.

An IR trial with center=.25,width=.08 gives R≈4.86×10⁻8 GeV, **below the previous 11-global-polynomial Ritz value 8.56×10⁻8 GeV**. This confirms the global polynomial basis was insufficient to characterize all physically relevant low-momentum disturbances. It also highlights how a positive self-adjoint collision form can have slower localized trial perturbations than relatively smooth thermal energy perturbations.

**Scientific scope**: masses, rates and T here are a demonstrative broken-EW χe contact EFT benchmark; the bath is fixed FD, no neutrino/UV mediator/comoving expansion/gstar(T) or stimulated massive Bose channel. No connection to physical 86-eV DM relic or CLASS is claimed.

[Code](localized_pauli_packets.py), [tests](../tests/test_localized_pauli_packets.py), [measured JSON](results/pauli-localized-packet-2026-11-21.json). Tests verify naive quadrature failure, IR composite n and cutoff refinement, physical G² scaling, positive finite susceptibility and parameter rejection.

**Next independent goal**: push localized packets to smaller width and radial p≈0 with a rigorously convergent collision quadrature (adaptive angular integrations and endpoint/cutoff estimates); if a spectral gap is conjectured absent, prove it from analytic bounds rather than finite Rayleigh trials. In parallel map actual χe collision event rates into a conservative Pauli+finite-SM kinetic transition network, not the toy arbitrary Γ.
