# Positive exponential filtering of closed-bath FD slope modes

The finite-bath entropy Hessian includes both visible and dark responses at fixed total energy. Compute this Hessian at the equilibrium associated with the actual scalar energy constraint. Normalize the positive finite-state timelike secant mobility by the equilibrium Hessian, replace each positive rate k by `2*tanh(dt*k/2)/dt`, and transform back. The conservative finite step then has the same total-energy constraint and chord rule for a modified total entropy production. No finite-step claim equates that modified production with the physical unfiltered production.

In two thermodynamic near-equilibrium tests, 0.30 and 0.68825 GeV, the measured fast-mode residual after a step with k_fast*dt=80 was below 2.49e-8 of the initial mode amplitude, while the exact linear coefficient is about 1.80e-35. The residual is measured from small nonlinear perturbations, not a resolved estimate of that tiny linear coefficient. Independent phase-space meshes agreed on mode amplitudes to better than 9.80e-11, and scaled total entropy balance defects in the sampled steps were at or below 1.05e-15.

Use the same EOS input as the [unfiltered tabulated bench](FD_PORTAL_TABULATED_BATH_ENTROPY_AND_MODES.md).
```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_tabulated_exponential_filter.py eos-range.dat
```

This is a fixed-energy numerical damping check, not a fully coupled expanding trajectory or particle-abundance result.
