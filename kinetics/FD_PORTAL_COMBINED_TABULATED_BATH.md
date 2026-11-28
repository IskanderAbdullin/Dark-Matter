# Elastic heat and timelike FD collisions in the same closed tabulated bath

Keep exactly the preceding tabulated visible energy and entropy; add direct finite-width vacuum spacelike elastic scattering on charged leptons. The two collision channels have the same heat sign convention (positive from dark to visible). The elastic channel carries no particle source, and its entropy term is the elastic heat multiplied by the temperature force. Verify the local entropy pairing and rate direction with two independent phase-space rules, without an interpolated scattering grid.

In three finite-displacement checks the coupled discrete-gradient update had error-halving ratios 4.0016–4.1119 against the differential system. At the T=0.3 GeV start, chosen near the timelike heat zero, adding elastic heat increased the net heat magnitude by a factor of 202 and contributed about 46.0% of local entropy production, while leaving the particle source at that exact state unchanged. The initial conditions in this comparison are controlled, not inferred cosmological initial conditions.

A separate equilibrium calculation across 13 temperatures adds the independently computed elastic heat coefficient as a positive rank-one mobility. At T=0.10 GeV the coefficient is only 16.23% of the timelike heat entry, but the slow rate grows by a factor 2.367; the fast rate grows by a factor 1.1215. At T=0.68825 GeV the factors are 1.01986 and 1.000247. Thus heat-only scattering can matter disproportionately for the slow coupled mode at moderate temperatures, even though it contributes no number source.

Use the authors' EOS input and normalization described in [the preceding tabulated-bath benchmark](FD_PORTAL_TABULATED_BATH_ENTROPY_AND_MODES.md); never redistribute the table.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error fd_portal_combined_tabulated_bath_entropy.py eos-range.dat
python -W error fd_portal_combined_tabulated_modes.py eos-range.dat
```

These are local closed-system checks with vacuum charged leptons only. No expanding history, complete interaction network, or relic density is implied.
