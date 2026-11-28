# Finite reservoir energy feedback in the local FD collision problem

Retain two independently evolving bath temperatures by making dark number and dark energy the step coordinates. Total energy is fixed; the visible energy is the remainder, with the explicitly ideal proxy `rho_vis = (10.75*pi^2/30) T_vis^4`. The negative total entropy at fixed energy has forces `(eta_D, 1/T_vis - 1/T_D)`. Recover an actual isoenergy state at the force averaged along the moment chord, and evaluate the full FD collision rate there. The conservative step evolves both temperatures and pairs that rate with the exact chord entropy identity. A separate 12-point chord rule, two phase-space orders, independent DOP853 reference and direct finite-difference entropy-gradient tests check the implementation. The tested error-halving ratios were 4.0019 to 4.1484. Both temperatures change, and the total energy of all saved test endpoints closes to roundoff.

At common temperature and chemical equilibrium, finite visible heat capacity adds the positive rank-one term `1/(T_vis^2 C_vis)` to the energy-energy entry of the negative entropy Hessian. A separate executable checks the isoenergy force derivatives by direct perturbations and resolves both collision modes in two independent quadratures. At fixed temperature and ideal visible effective degrees of freedom 10.75, the fast rate is 1.3447 to 1.4086 times the fixed-reservoir value across the four sampled temperatures, while the slower rate increases by 0.7% to 5.5%.

The ideal reservoir is a controlled thermodynamic model, not the tabulated Standard Model equation of state. These are closed, homogeneous collision-only diagnostics, not expanding-universe, decoupling, or abundance calculations.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest -q test_fd_portal_finite_bath_discrete_entropy
python -W error fd_portal_finite_bath_discrete_entropy.py
python -W error fd_portal_finite_capacity_modes_reference.py
```
