# Full finite χ/e Pauli energy-lattice dynamics: bounded all-quartet test, not physical cosmology

Code: kinetics/paired_quartet_network.py ; tests: tests/test_paired_quartet_network.py ; reproducible runner: python kinetics/verify_paired_quartet_network.py --output outputs/quartet-time.json ; compact results: kinetics/results/paired-quartet-network-2026-11-21.json.

## Scope and falsifiable experiment

This stage moves beyond the single-quartet update of PAULI-CHI-E-QUARTET-AUDIT.md by evolving **all** allowed discrete χ+e ↔ χ′+e′ energy-conserving quartets in a finite lattice. Crucial qualifier: the positive W tensor is computed from the actual tree vector-contact full-Pauli 5D microphysical scattering, BUT uses a declared arithmetic pair-average (Wraw+Wreverse)/2 to repair raw finite-grid reversibility. **This is a modified equilibrium-consistent lattice model**, NOT the exact original continuum nonlinear χe collision integral or a physical γ_el for relativistic BEST. Every observable retains this qualification.

For every unordered energy-matched quartet i,a↔j,b, we use the independently unit-tested fully implicit Pauli transfer:
  ∆N=∆t Wpair[F_new−R_new], with
  F=fi ga(1−fj)(1−gb) and R=fj gb(1−fi)(1−ga).
The scalar root is constrained by the exact four capacity bounds, so 0≤all f≤1 without a poststep clip. It conserves χ number and electron number separately, and χ+electron total energy exactly because Eχ_i+Ee_a=Eχ_j+Ee_b on the common lattice. The local Fermi entropy increases. The network updates ALL positive-energy-transfer unique quartets by sequential Lie splitting, with forward/reverse ordering checked independently.

No external fixed-T electron reservoir exists in this test: χ and electron modes **both evolve** and exchange energy self-consistently within the discretized finite phase-space. That does NOT include QED plasma thermal masses, other SM fields, lepton chemical processing, Friedmann expansion, physical χ1↔χ2+φ Bose enhancement or electroweak UV completion.

## Measured time evolution, independent collision benchmark

Fixed-volume test: Mχ=0.1 GeV, G=0.001 GeV^-2, Tχ,initial=1.2 GeV, Te,initial=1.0 GeV, μχ=μe=0, 14 χ energy nodes, 17 electron energy nodes, spacing δE=0.5 GeV, three collision angles with 8 nodes each. At this small cutoff, this is a **numerical algorithm benchmark**, not a realistic physical long-time plasma. The source has 1,092 positive-rate, energy-changing quartets.

The pre-step instantaneous full-Pauli **discrete** heat source evaluated by a completely separate evaluation of the tensor is Qχ(initial)=−4.96747419×10^-8 GeV^5. For dt=1000 GeV^-1, the time-integrated result ∆Eχ/dt=−4.96575086×10^-8 GeV^5, a −0.0347% relative deviation; dt=10000 gives about −0.346% and dt=100000 about −3.36%. This is an actual timestep truncation effect, not a violation of energy conservation.

For a much longer total interval dt=2×10^6 GeV^-1, **forward-split** χ energy changes are:
- 1 substep: −0.05791998857 GeV^4; entropy gain +0.00619333136 GeV^3.
- 4 substeps: −0.05571543344 GeV^4; entropy gain +0.00604842272.
- 16 substeps: −0.05565715219 GeV^4; entropy gain +0.00603205902.
- 64 substeps: −0.05567090746 GeV^4; entropy gain +0.00602906076.

All these runs preserve separate χ and e number to ~10^-16 GeV^3, total energy to ~10^-15 GeV^4, and Pauli occupation bounds. But **this does not prove time convergence**: the ordering difference between forward and reverse splitting is 1.44282×10^-4 GeV^4 (one substep), 2.26681×10^-5 (4), 5.43156×10^-6 (16), and 1.35004×10^-6 (64), showing first-order-like reduction but residual error at every finite dt. A thermal common-T FD state with independent μχ=+0.15, μe=−0.2 remains stationary to <10^-12 for dt=10^7, as required by local detailed balance of the MODIFIED Wpair lattice.

**Negative lesson:** unit tests of exact N, E and positive entropy alone CANNOT certify physical accuracy: both the original 5D→finite-bin reverse-rate defect and temporal split error remain quantifiably nonzero. Retain as falsification controls; do not tune an ad hoc global BGK γ to hide them.

## Next decisive physics

Compute physically identical forward/reverse on-shell scattering phase-space measures or bounded within-cell detailed-balance approximation without simply averaging two rate tables. Next refactor the full operator into an efficient sparse χ/e collision network with measured energy-bin, angular, time, cutoff and stochastic errors and compare its nonthermal heat/history with the continuous 5D weak/strong χe source. A true gauge-consistent χ1–SM portal, full e± and neutrino thermal bath, physical 1↔2 Bose daughter kinetics, gstar(T), FRW freezeout, qmin/zero-mode and observational CLASS comparison remain separate unsolved milestones. No physical 86-eV relic, BEC, ΩDM or CLASS-ready is claimed.

Literature context: Binder et al. (2017) https://arxiv.org/abs/1706.07433 and Binder et al. (2021) https://doi.org/10.1140/epjc/s10052-021-09357-5 on momentum-resolved Boltzmann versus reduced kinetic closures. Neither supplies an independent numerical validation of this finite-bin scheme.
