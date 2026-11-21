# Microphysical 5D Pauli χe scattering to finite bins: a negative reversibility test

Status: source and independently run sandbox numerical experiment. The GitHub CI must be checked separately. This is an electron-only broken-electroweak vector-contact EFT with Mχ=0.1 GeV, T=1 GeV, G=0.001 GeV^-2, gχ=ge=4, zero chemical potentials and massless e±. This DOES NOT solve an arbitrary-f Boltzmann equation, gauge/UV plasma, finite-temperature Γ(T), 86-eV Bose relic, BEC, ΩDM, or CLASS.

## Pre-registered hypothesis and independent reference

Finite-volume energy-bin projection of the exact positive 5D χe scattering events should conserve the event's FIRST χ energy moment, but it need not conserve local FD detailed balance or the quadratic Dirichlet collision form. The independent reference is the earlier full continuous 5D full-Pauli FD energy-susceptibility Dirichlet D11=3.02789980222e-7 GeV^4 from LINEARIZED-PAULI-EIGENMODES.md, computed by a separate Gauss-p formulation. A pairwise reversible closure is NOT assumed to be microscopically exact.

The implemented microphysical event weight for an initial χ radial bin i is

    dL = V_i [g_e/(2π)] k²dk dμ v_M (dσ/dΩ_CM) du_CM <dφ_CM/2π>
         × fχ(E_i) fe(k) [1-fχ(E′)] [1-fe(k′)]

where V_i=gχ [p_upper³-p_lower³]/(6π²) is the EXACT bin-state capacity, v_M=1-p μ/E, dσ/dΩ_CM=|M|²/(64π²s), and the outgoing E′, k′ come from the exact 2-body Lorentz boost. Thus electron recoil and all FOUR equilibrium Pauli factors are included. Kinematic tails above the last χ node are counted as omitted weight, not silently projected into the grid; at pmax=kmax=22 GeV it is about 5.2×10^-7 of the incoming event weight in the demonstrated thermal benchmark.

Interpolate E′ between neighboring energy nodes E_j,E_(j+1) with positive barycentric weights B_j(E′). This conserves the event's FIRST outgoing energy exactly when E′ lies within the grid, but it does not conserve its squared energy displacement or the equilibrium reverse rate.

## Independent five-dimensional numerical measurements and failures

For χ energy nodes formed by p_i=22[i/(N-1)]^1.6 (GeV), finite radial p³ cell volumes, and 12 Gauss k nodes, 10 incoming angle, 10 CM polar and 10 CM azimuth nodes:

| Number of χ nodes | Physical unprojected D, GeV^4 | Projected D, GeV^4 | Excess projected quadratic dissipation |
|---:|---:|---:|---:|
| 16 | 3.12036008×10^-7 | 3.34439139×10^-7 | +7.1813% |
| 24 | 3.06692466×10^-7 | 3.16088725×10^-7 | +3.0653% |

For 64 nodes and 16 nodes in EACH of the four other dimensions, physical D=3.03301126×10^-7 GeV^4, projected D=3.04545928×10^-7, a +0.4119% interpolation excess. At χ=128 with 20 electron energy and 16 angular nodes, physical D=3.02906773×10^-7 GeV^4, within +0.0386% of the independent full-continuum n24/cut32 calculation. The projected excess is then +0.1013%. These convergence numbers test ONE thermal equilibrium energy-mode functional, not full arbitrary-f spectra, kinetic eigenmode gap or cosmological physics.

**Exact mathematical source of bias:** because ∑_j B_j=1 and ∑_j B_j E_j=E′, the cell-projected quadratic loss satisfies

    ∑_j B_j(E′)(E_j-E_i)^2
      = (E′-E_i)^2 + ∑_j B_j(E′)(E_j-E′)^2.

The final term is a nonnegative Jensen variance. Therefore merely conserving each projected event's average energy does NOT ensure that the linearized scattering relaxation spectrum is correct. All reported positive excess values are genuine finite-grid discretization errors, not extra physical scattering.

**Second negative result — raw microscopic projected rates are not paired.** Define L_ij by the positive physical 5D event weight times B_j(E′). At χ=64 nodes with 16-point other dimensions, ∑|L_ij-L_ji|/(2∑L_ij)=1.5361%; the raw FD stationary number-source L1 defect normalized by ∑L is 1.3245%. At 24 angular nodes the unpaired asymmetry decreases to 0.5591%, despite an almost unchanged physical D. Even at exact common T the raw finite-bin projection produces a spurious net χ heat +1.69647×10^-10 GeV^5 for χ64/angular16, whereas the exact continuum FD collision is identically zero. This bias is NOT removed by angular convergence of the energy quadratic alone.

The deliberately declared approximation L_rev=(L_raw+L_raw^T)/2 exactly restores the finite-bin FD stationary point, but CHANGES the physical transition graph. Its reversible pair coefficient Γ_ij=L_rev,ij/[f_eq,j(1-f_eq,i)] is derived from the physically integrated equilibrium events plus that modification. It is NOT an independently derived full nonlinear rate at an arbitrary bath temperature. The centered numerical ε=10^-4 logit energy perturbation gives graph dQ/dε=-3.0454597360e-7 GeV^4, agreeing with minus its projected D to 1.50×10^-7 relative; this is a self-consistency test of the declared reversible approximation, not proof of its physical accuracy.

The prior implicit finite-bath pair integrator (fermion_finite_bath_pairs.py) can consume one such Γ calibrated at T0=1 GeV and preserve Pauli, charge, energy and nondecreasing total entropy. However, keeping Γ fixed when T_SM changes is not a recomputation of the full electron FD Pauli scattering and cannot validate actual evolving plasma kinetics. The evidence is therefore a constructive physical 5D equilibrium collision projection AND a quantified failure to satisfy required finite-bin physical detailed balance without an uncontrolled correction.

## Reproduction and literature

Code: kinetics/microphysical_equilibrium_graph.py
Unit tests: tests/test_microphysical_equilibrium_graph.py
Compact reference: kinetics/results/microphysical-bin-graph-2026-11-21.json
Runner: python kinetics/verify_microphysical_bin_graph.py --output outputs/real-5d-bin-audit.json

Relevant external approaches: Binder, Bringmann, Gustafsson & Hryczuk (2017), Early kinetic decoupling of dark matter, https://arxiv.org/abs/1706.07433 ; Binder et al. (2021), Dark matter relic abundance beyond kinetic equilibrium, https://doi.org/10.1140/epjc/s10052-021-09357-5 . They motivate phase-space Boltzmann rather than assuming a universal rate; neither is used as independent numerical validation of our new graph.

## Next necessary physics

Create an *energy-resolved paired χ+e transition quadrature* (joint collision quartet), with exactly matching initial/final total energies and a common positive microphysical rate on each pair. Either use a joint χ/e fine energy lattice with exact energy-conserving reaction quartets and full Pauli, or derive a controlled within-cell equilibrium closure with independently bounded bias. Verify the reverse rate without post-hoc symmetrization, preserve individual χ/e charge and total χ+e four-momentum, and only then evolve with finite SM entropy and the pinned BEST 1↔2 Bose channels. Do not promote Γ at one reference T to an all-temperature γ_el or cosmological result.
