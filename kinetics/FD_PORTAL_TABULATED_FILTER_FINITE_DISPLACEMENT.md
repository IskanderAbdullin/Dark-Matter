# Finite-displacement refinement of the closed tabulated FD filter

Use the same three fixed-energy, tabulated-visible-bath starts and compare the filtered update against an independent DOP853 solution of the **unfiltered** timelike moment ODE. The finite-state positive filter modifies the discrete collision mobility but does not modify this continuous reference. Four-step endpoints agree on the independent phase-space grids. The error-halving ratios over 2/4/8 steps are 4.0045–4.0812.

At eight steps, normalized moment errors at initial visible T=0.05, 0.30 and 0.68825 GeV are 8.88e-7, 3.11e-5 and 1.09e-5 with the filter. The corresponding unfiltered values are 5.98e-6, 2.58e-5 and 4.60e-5. The filter improves two of these sampled outcomes but is not uniformly more accurate. Its total-entropy identity is the modified discrete identity, not equality to the physical finite-step entropy source.

Use the same local authors' EOS normalization as in [the tabulated FD benchmark](FD_PORTAL_TABULATED_BATH_ENTROPY_AND_MODES.md).

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_tabulated_filter_finite_displacement.py eos-range.dat
```

These tests use closed homogeneous timelike collision substeps, not an expanding thermal history or a relic density estimate.
