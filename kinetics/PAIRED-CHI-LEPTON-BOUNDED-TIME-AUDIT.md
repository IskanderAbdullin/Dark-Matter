# Nonlinear Pauli-bounded χ–lepton 2→2 time pilot on finite phase-space energy cells (2026-11-22)

## Question and acceptance criterion set before the experiment

The prior on-shell finite-vector 5D χℓ event weak heat is now projected into separate energy bins using matching Fermi-logit interpolation and barycentric deposition. That achieves exact instantaneous species number, exchanged energy and quantum H-theorem, but the projected source had NOT been evolved as a true changing fermion population. The previously disproved naive interpolation had negative binned entropy for an admissible oscillatory Fermi distribution.

Test a genuine coupled **nodewise nonlinear time pilot** recomputing the full 5D Pauli source at every step, without FD thermal resets. Criteria: each node must remain 0<f<1, all accepted steps increase actual Fermi entropy, χ and lepton numbers and total energy remain constant in a closed 2→2 system. An adaptive full-spectrum step-doubling check must reject excessive steps and respond to tighter tolerances. Report the cost, phase-space measure, time-unit limitations and the still-uncontrolled momentum-grid error.

## Physical radial cells, source and algorithm

For each species with rest mass m and historical numerical spin factor g=4, use energy nodes E_i and the clipped physical radial cell [L_i,U_i] within [m,Emax]. Define the exact **zeroth density-of-states measure**

    w_i = g/(6π²)*[(U_i²−m²)^(3/2)−(L_i²−m²)^(3/2)] > 0.

These cells exactly partition the physical volume g pmax³/(6π²). The first energy moment remains collocated at the bin node; this is NOT an exact continuum cell-energy or scattering solver. Both g factors here use the historical 4×4 5D normalization, not correct physical multiplicities for every SM neutrino flavor.

At every time step, reconstruct the event occupations from the CURRENT Fermi node logits, with the SAME nonnegative barycentric weights as signed event deposition, and recompute the relativistic finite mediator 2→2 Pauli forward-minus-inverse integral. The two nonlinear node equations are

    d fχ_i / dt = (d Nχ_i / dt) / wχ_i,
    d fℓ_j / dt = (d Nℓ_j / dt) / wℓ_j.

These inherit source-level conservation Nχ, Nℓ, energy Qχ+Qℓ=0 and exact instantaneous bin H theorem Sdot=Σ W(F−R)ln(F/R)≥0. For each candidate dt, try one full explicit Euler step AND two independently reevaluated half-Euler steps. Accept the TWO half-step solution only if all intermediate and final nodes stay strictly within (0,1), actual finite Fermi entropy does not decrease, and the maximum over both species of |one-full − two-half| is at most atol. Otherwise reject and halve dt. The next step is estimated from the local error with guarded growth. There is no artificial clipping or post-step charge repair.

**Important:** This numerical time method remains FIRST ORDER despite the step-doubling error control; the accepted local discrepancy does not certify a global error bound or physical kinetic rate.

## Reproduced compact numerical experiment

Fixed illustrative dimensionless setting Mχ=.1, mℓ=0, G=.001, mV=.1, Te=1, kinetic cutoff12, incoming nrad=nmu=4, outgoing polar/azimuth=12×12, only TEN energy intervals per species. The total numerical time is 10^10 in units defined by the arbitrarily normalized 5D rate, initial proposed step 10^10. This is not cosmological time or physical Γelastic/H.

Two distinct synthetic initial distributions, subsequently evolved without any thermal reset:

- Oscillatory fχ=expit(−E+0.3 sin 4E), fℓ=expit(−ω+0.3 cos 4ω); the old H-violation control.
- Smooth fχ=expit(−E/.9+0.15 exp(−E/.9)), fℓ=expit(−ω+0.05 exp(−ω)).

| Initial shape | atol | Accepted / rejected | Maximum local spectral step-doubling error | Total entropy gain |
|---|---:|---:|---:|---:|
| Oscillatory | 10⁻⁶ | 4 / 2 | 9.70×10⁻⁷ | +4.21465×10⁻⁴ |
| Oscillatory | 10⁻⁷ | 11 / 7 | 9.77×10⁻⁸ | +4.21146×10⁻⁴ |
| Oscillatory | 10⁻⁸ | 33 / 15 | 9.79×10⁻⁹ | +4.21049×10⁻⁴ |
| Smooth | 10⁻⁶ | 2 / 1 | 8.80×10⁻⁷ | +9.94502×10⁻⁵ |
| Smooth | 10⁻⁷ | 7 / 4 | 9.75×10⁻⁸ | +9.92979×10⁻⁵ |
| Smooth | 10⁻⁸ | 24 / 12 | 9.79×10⁻⁹ | +9.92540×10⁻⁵ |

The independent sandbox reconstruction measured weighted χ number, lepton number and total energy residuals no greater than 4.5×10⁻¹⁶, with strictly positive entropy gain per accepted substep and all Fermi nodes within physical bounds. Tighter tolerances cost more real 5D collision evaluations and visibly refine the first χ node (oscillatory final first occupation approximately .503052630, .503053576, .503053863).

This establishes a useful *finite-grid nonlinear time pilot* and exercises real adaptive rejection. It does not establish that the coarse-grid evolution is an accurate physical scattering rate. Prior independently measured oscillary spectral entropy alias −93.4% at 10 energy intervals survives untouched. Furthermore nr4, outgoing 12×12 and UV cutoff12 have not been converged.

## Physical caveats and independent context

True SM leptons require correct spin/flavor and ν vs e± gauge coupling accounting. Both χ1 and χ2 must be evolved with proper abundance under the simple real-neutral-φ Yukawa equal-charge U1′ branch, along with the on-shell mediator, finite SM ν/e±γ bath and Bose/Pauli Yukawa χ1↔χ2+φ. None of these processes or FRW H have been coupled here. No physical Γelastic, frozen 86-eV f(q), ΩDM, macroscopic BEC or CLASS-ready result is established.

Methodological context, not imported code: [entropy-stable nonlinear Boltzmann, JCP 463 (2022) 111289](https://www.sciencedirect.com/science/article/pii/S0021999122003515) and [adaptive far-from-equilibrium Boltzmann numerical scheme, CPC (2021)](https://www.sciencedirect.com/science/article/abs/pii/S0010465521000254).

Reproduce from the permanently archived [source](paired_chi_lepton_bounded_time.py), [tests](../tests/test_paired_chi_lepton_bounded_time.py), [JSON](results/paired-chi-lepton-bounded-time-2026-11-22.json), [dated stage](../archive/stages/2026/11/2026-11-22-chi-lepton-bounded-time.md):

    python -m pytest -q tests/test_paired_chi_lepton_bounded_time.py
    python -m kinetics.paired_chi_lepton_bounded_time

Next P0: verify the real GitHub Scientific CI, independently scan at least three energy, incoming, outgoing and cutoff resolutions, error-budget the full chi/lepton spectra and true 5D weak/strong equality, then couple χ1+χ2+φ Bose and complete SM/Z′/FRW physical parameters. Never promote this first-order ten-bin pilot into a terminal dark matter claim.
