# Time-ordered response of the late expanding FD history

In the full comoving variables (dark number, dark energy, total energy), split the derivatives with respect to log(a/a_start) into collision and expansion contributions. Both receive the physical temperature and chemical potential reconstructed from the densities; H varies with total energy. The collision block for perturbations in dark number and energy at fixed total energy reproduces the independent fixed-energy implementation within 3.46e-14. The full local matrices at the twelve stored timelike/combined points have three negative real eigenvalues; that alone does not imply contraction of a nonautonomous tangent.

For the combined solution between x=.96 and 1.20, integrate the time-dependent matrix in one fixed initial comoving basis. The largest singular value of the endpoint response is **0.9871126** when both ends are expressed in the starting background units. Rescaling the output by the final background densities gives **1.0187359** when the response is measured in endpoint fractions. Direct perturbed trajectories give 1.0187359 in that same endpoint norm. A unit direction for the latter response is (+0.439777, +0.260076, +0.859626) in initial fractional number, dark energy and total energy. If the initial total comoving energy is fixed, restrict the input to the first two columns while keeping all three output fractions. The combined gain drops to 0.9767519 (direct perturbations: 0.9767519); thus the full response above unity requires allowing the initial total energy to vary. These gains are properties of the chosen observable normalization over one interval, not physical stability or instability statements.

In the independent timelike calculation the corresponding gains are 0.9876592 and 1.0165415. The largest matrix-element difference from the combined case in their respective initial background normalizations is 1.37e-02. The timelike ordered response agrees with its own direct perturbations within 2.47e-07. These comparisons change the reaction channels, not the interpretation of the basis.

The 13/25/49 sampling comparisons have maximum elementwise differences 1.10e-06 and 7.60e-07. Independently integrating the six nonlinear perturbations reproduces the ordered response within 3.03e-07. The log determinant matches the integrated trace within 5.55e-16. The test suite also checks an exactly soluble noncommuting example and verifies that simple collisionless radiation has no response gain when measured in each endpoint's fractions.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fd-portal; e=eos-range.dat
python -W error -m unittest -v test_fd_portal_direct_comoving_variational_propagator
python -W error fd_portal_direct_comoving_variational_jacobian.py "$e"  "$p-direct-expansion-055-096-$d.json" "$p-direct-expansion-096-120-$d.json"
python -W error fd_portal_direct_comoving_variational_propagator.py "$e"  "$p-direct-expansion-096-120-$d.json" "$p-direct-expansion-055-096-$d.json"
```

Repeat the second integration with `--operator timelike` to reproduce the timelike receipt. Use the same normalized public EOS tabulation. No new hadronic channel, general freeze-out criterion, abundance constraint or transfer function is implied by this numerical test.
