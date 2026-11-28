# Positive finite-step mobility filter for a local FD collision step

Scale number and energy by the constants `S=diag(2,T_SM)`, so the dark force vector is conjugate to `q=S^{-1}(n,rho)`. Its equilibrium susceptibility `C=dq/du` is symmetric positive definite. At the force averaged along a moment chord, the exact FD number/heat flux has secant mobility `H`, giving the symmetric normalized matrix `K=C^{-1/2} H C^{-1/2}`.

For a positive step replace the eigenvalues `k` of `K` with `2*tanh(dt*k/2)/dt`, and transform back by congruence with `C^{1/2}`. Every filtered eigenvalue is positive, so the chord chain rule still gives a nonnegative entropy production for the modified flux. In the equilibrium linearization the resulting Cayley factor is exactly `exp(-dt*k)`, avoiding the stiff sign reversal. The map tends to the unfiltered mobility as `dt` tends to zero, with an O(dt^2) perturbation at fixed `k`.

The filtered finite-step entropy is **not** the original physical entropy rate. The test records their ratio as well as the balance residual and verifies damping with positive and negative finite perturbations on two independent phase-space meshes. This is a local numerical candidate, not a full expanding history or a general stiffness guarantee.

On the fast mode at z=80, the largest absolute normalized nonlinear response was 2.81e-8 across the three temperatures at the smaller perturbation. The filtered step entropy rate was close to 2.5% of the original physical rate at the same averaged force. Discrete entropy identity error stayed below 1.59e-15. The largest between-mesh amplitude difference across all tests was 1.02e-10.

Run from `kinetics`:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_reservoir_exponential_filter_reference.py
```
