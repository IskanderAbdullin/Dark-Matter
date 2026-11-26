# Pauli χ+e physical scattering on a common energy lattice: microreversible QUARTETS

Status: independently tested 5D physical tree-level χe matrix element mapped to a **finite discrete** common energy lattice with a declared reversible pair closure. It is not the unmodified continuum nonlinear collision kernel, not a full multi-quartet time integration or a cosmological dark matter solution. Source: kinetics/paired_energy_quartets.py; tests: tests/test_paired_energy_quartets.py; compact reproducible measurements: kinetics/results/paired-energy-quartets-2026-11-21.json.

## Pre-registered hypothesis and physics

A common energy lattice with χ_i energy Eχ_i=M+(i+1/2)δ and massless e_a energy k_a=(a+1/2)δ allows an individual discrete quartet (i,a) ↔ (j,b), b=i+a−j, to satisfy exact center total-energy equality Eχ_i+k_a=Eχ_j+k_b. Each species has exact radial number capacity Vχ_i=gχ/[6π²](p_{i+1}³−p_i³), Ve_a=ge/[6π²](k_{a+1}³−k_a³). This is an **energy-centered finite-bin approximation**; continuous kinematics is evaluated using true on-shell χ and electron incoming momenta, exact CM outgoing boost, vector-contact spin-averaged |M|², Møller flux and all the physical outgoing energies.

An incoming event with continuous final Eχ' is deposited into the two nearest χ energy centers with positive linear barycentric fractions, while the electron destination b=i+a−j enforces exact discrete total energy. If j or b falls outside the domain, its weighted event rate is discarded and separately measured. The bare physical projected tensor Wraw_(ia,jb) is positive and has units GeV^4, without any χ/e occupancy factors; physical occupancies are inserted in BOTH Pauli directions when computing rates.

**Negative microscopic reversibility control:** this positive 5D event projection does NOT yield Wraw_(ia,jb)=Wraw_(jb,ia) at finite grids. Define Wpair as their average and label this an explicit **model-changing symmetric rate closure**, not "the true microphysics". With that declared closure, individual pair source is

    J_(ia,jb) = Wpair_(ia,jb) [
      fχ_i fe_a (1-fχ_j)(1-fe_b) −
      fχ_j fe_b (1-fχ_i)(1-fe_a) ].

For each energy-matched quartet it exactly conserves χ number, electron number and χ+e energy and produces nonnegative entropy rate J ln(F/R). At any common T, independent χ and electron chemical potentials give F=R pointwise; no external bath assumption is required for this discrete identity. At f=0 or 1, the diagnostic entropy integral is withheld until its limiting value is handled (implicit single-quartet update itself accepts fermionic bounds).

## Measured independent collision/heat comparison

Demonstration model: Mχ=0.1 GeV, Te=1 GeV, G=0.001 GeV^-2, gχ=ge=4, zero μ, electron mass neglected, reference bath Te=1. For χ/e grid δ=0.5 GeV, Nχ=36, Ne=54 and nμ=nCM=nφ=10 (χ upper center 17.85 GeV; electron up to 26.75 GeV), the physically integrated one-time Pauli source, with the explicitly paired rate closure, yields:

- χ colder, Tχ=0.9 GeV: Qχ=+2.01553821909×10^-8 GeV^5, versus separately integrated continuous 5D Gauss weak reference +1.98528093330×10^-8 GeV^5, relative difference +1.52%.
- χ hotter, Tχ=1.1 GeV: Qχ=−4.48214575005×10^-8 GeV^5, versus separate continuous 5D Gauss weak reference −4.43252672179×10^-8 GeV^5, relative difference +1.12% in magnitude.
- Qχ+Qe=0, χ and electron number sources each integrate to ≈0, and positive entropy for unequal temperatures; at common T, arbitrary independent μχ,μe, the paired source vanishes at floating-point precision.
- The raw (unpaired) finite lattice **still** has a nonzero false equilibrium Qχ=+5.79506×10^-12 GeV^5 and 3.6935% Pauli-weighted forward/reverse tensor asymmetry at angular n10. This is a real numerical failure of the finite-bin projection, masked by using the rate-averaged closure.
- About 4.16×10^-5 of FD-weighted event rate leaves the supported quartet lattice. A coarse χ24/e48, δ=.5 grid has ~0.187% FD-weighted omission and hot χ Q=−4.2203×10^-8 GeV^5. **Grid/domain convergence remains incomplete**; the 1–2% agreement is a benchmark at the specified resolution, NOT a universal accuracy bound.

A smaller χ30/e30, δ=.4 lattice has Pauli-weighted raw tensor asymmetry 2.8954% for angular n12, reduced to 1.5264% for n16 at approximately fixed projected Dirichlet energy moment. Like the previous raw χ-only graph, this demonstrates that a converged heat moment alone does not guarantee local reversibility.

## Real-time two-species microphysical-rate SUBSTEP (not full time evolution)

The single-quartet implicit collision step uses the positive Wpair from the actual 5D scatter evaluation and solves for a scalar transferred NUMBER ΔN. It updates all four fermionic occupations together through their physical capacities:

    fχ_i(new)=fχ_i(old)−ΔN/Vχ_i
    fe_a(new)=fe_a(old)−ΔN/Ve_a
    fχ_j(new)=fχ_j(old)+ΔN/Vχ_j
    fe_b(new)=fe_b(old)+ΔN/Ve_b
    ΔN=dt Wpair (F(new)−R(new)).

The exact admissible ΔN bracket enforces 0≤all four f≤1 for arbitrarily stiff dt. For a sampled nonstationary quartet at dt=10^8 GeV^-1, the measured χ/e charge residuals and total energy residual are ≪10^-12 in their native units; total fermionic entropy increases. No post-update clipping, hidden χ thermal reset or external finite heat reservoir is used: electron occupations evolve on this quartet.

**But:** only ONE physically sampled quartet is advanced at a time. No whole coupled χ/e Boltzmann history, simultaneous collision network convergence, FRW expansion, physical 1↔2 quantum/Bose emission, realistic SM e/ν/muon bath or UV-complete mediator has been solved. A long-time 86-eV freezeout, BEC, ΩDM, true kinetic γ_el and CLASS predictions are NOT justified by this step.

## Next independent scientific gate

Derive **directly paired** event/phase-space quadrature on the identical initial/final collision manifold without global a posteriori W symmetrization. Measure how the reverse asymmetry, grid and angular errors respond as the lattice δ→0 at fixed sufficiently broad p/k support. Derive an efficient sparse microreversible χ/e quartet operator and a provably Pauli-bounded global splitting or implicit integrator. Compare its time-dependent Qχ and entropy, fermion charges and heat spectra with the separate continuous full-Pauli 5D source under q, angular, energy-cutoff and dt refinement.

Only afterwards consider electroweak lepton-doublet/ν channels and UV completion and coupling to pinned BEST χ1↔χ2+86-eV φ and realistic gstar/Friedmann. Related independent kinetic-method literature: Binder et al. (2017) https://arxiv.org/abs/1706.07433 ; Binder et al. (2021) https://doi.org/10.1140/epjc/s10052-021-09357-5 ; these are NOT claimed as independent validation of this finite-bin approximation.
