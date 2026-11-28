# Local inverse of the actual bath-dependent FD number/heat response

Hold the visible temperature and the finite timelike lepton portal fixed.
Generate target rates at five prescribed force points near the accepted
nonequilibrium samples. Invert with the analytic **physical** rate
Jacobian, including the change of event activity with both forces.
Redo the inverse on a different positive quadrature and with a
central-difference Jacobian from the independent rate evaluator.
The analytic inverse recovers every prescribed target in three accepted
steps on both meshes, with normalized residual below 1.1e−15.
The finite-difference solutions agree in the same force neighborhood.

For comparison, invert the very same targets using the convex
frozen-activity potential of the starting state. Although that
potential can reproduce its own frozen targets, the resulting
forces miss the actual required rates. The normalized physical
residual ranges from 0.0149% to 1.152% here, using the diagonal
Jacobian scaling described in the code. Neither the secant nor
the frozen Hessian can generally substitute for the actual tangent.

The exercise demonstrates a local numerical inverse for these
prescribed controls, not a global inverse, a bath trajectory, or a
cosmological stability or abundance prediction.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_physical_local_inverse_reference.py
```
