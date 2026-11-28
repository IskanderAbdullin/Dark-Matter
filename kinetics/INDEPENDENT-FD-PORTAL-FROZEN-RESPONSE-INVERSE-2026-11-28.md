# Conditioned inverse of the positive frozen-activity FD response

Freeze the accepted microscopic forward/reverse activity at each
local visible temperature and pair of forces. The convex test-force
potential from the preceding control then has a symmetric positive
Hessian H. Recover both the original force and a prescribed nearby
test force from their generated flux pairs with damped Newton and
diagonal scaling. Confirm a vanishing primal-dual Fenchel gap at the
inverse target and a strictly positive one away from it.

In this fixed activity problem, keeping the number flux unchanged
gives a heat-only derivative H22−H12²/H11. Directly invert the two
oppositely displaced heat-only targets; the finite difference recovers
that estimate on both independently constructed quadrature meshes.
At visible T=.64135 GeV it is 0.8102% of the unconditioned heat
diagonal. At the sampled .30 GeV point it is 17.29%. Correlation with
number flux therefore matters to the inverse in these held-state
controls. The five prescribed nearby force pairs are recovered in
two to six Newton steps, with scaled residuals below 8.3e−14.

The inverse holds the event activity fixed. It is not an inverse of
the actual bath-dependent flux field with its nonzero curl, and does
not establish an evolving physical heat transfer or a stability bound.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_frozen_response_inverse_reference.py
```
