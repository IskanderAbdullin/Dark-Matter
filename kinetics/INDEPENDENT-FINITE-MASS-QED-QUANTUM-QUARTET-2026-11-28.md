# One reversible quantum event with both charged and Bose cells evolving

Replace the held charged coefficients of the [Bose graph control](INDEPENDENT-FINITE-MASS-QED-BOSE-EDGE-NETWORK-2026-11-28.md) by two bounded fermion occupations. This fixture isolates one allowed event P+γ ↔ b+s with a single fixed amplitude scale. It does not project an actual plasma collision integral.

One scalar event extent I changes the capacity-weighted inventories in the ratio
`(-1,+1,-1,+1)` for parent, daughter, photon and scalar. This exactly preserves both the charged count and the Bose-pair count. It preserves their combined energy when the selected energies satisfy E_P+ω=E_b+E_s. For the positive event the Pauli-Bose statistical rate is

```text
f_P (1−f_b) n_γ (1+n_s) − (1−f_P) f_b (1+n_γ) n_s.
```

Choose the initial direction from the separate nonnegative directional products. Interior states can equivalently use the log affinity; for an exactly empty or completely blocked cell, use its one finite direction without evaluating an undefined logarithm. The admissible extent stops at the first occupation barrier. Along that segment the log forward-to-inverse ratio strictly decreases, as each occupied-to-unoccupied log ratio changes in the opposing direction. This gives a unique balancing point within the segment and nonnegative total entropy production. If both products are exactly zero the isolated edge is stationary. If one is positive, occupations move into the open interval before the balancing point. The implementation finds the balance via log odds, then independently factors the finite-degree extent polynomial as (r−x)H(x), where x is the dimensionless extent and r is its equilibrium value. The variable −log(1−x/r) removes the logarithmic singularity from the time integral, making the positive time-of-flight equation easier to invert. Very late states within exp(−65) of the limit are returned at float64 equilibrium accuracy.

For 45 controls the time-of-flight result agrees with independent DOP853 to 2.09e−14 absolute, two partial steps agree to 8.89e−16, and the per-cell event ledger agrees to 6.13e−16. The energy ledger closes to 1.39e−17. Thermal Fermi/Bose occupations with one common temperature and conserved fermionic chemical potential are unchanged. Both directions of transfer are included. In two large-occupation stresses with scalar baselines near 1e16 and 1e15, a separate scalar-extent integration reproduces event extents +.05709 and −.18983 to 6.11e−16 absolute. A direct scalar snapshot difference returns 0 and −.25 instead. The scalar baseline is a stress parameter; the capacities and other occupations are fixtures, not dynamical archive predictions.

Six additional Pauli-Bose boundary controls cover a full fermion, an empty fermion, an empty photon, an empty scalar and a doubly blocked edge. They agree with DOP853 at better than 1e−14 absolute, with no division by a vanishing occupation and no seeded population where both products are zero.

Reproduce with `python -W error finite_mass_qed_isolated_quantum_quartet_reference.py` from `kinetics/` and compare the [numerical receipt](results/finite-mass-qed-isolated-quantum-quartet-2026-11-28.json). This test assumes one reversible rate element; it does not choose dressed kernels, fit occupancies or solve a network of multiple physical channels.
