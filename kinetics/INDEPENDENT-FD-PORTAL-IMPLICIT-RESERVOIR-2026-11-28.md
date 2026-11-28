# Entropy-respecting implicit collision updates against a fixed visible reservoir

Use the finite TIMELIKE charged-lepton FD model and the independent
dark number and energy state functions (massive FD and massless
scalar). Hold the reservoir temperature fixed. Each implicit step
uses the true Jacobian, including actual population dependence,
and solves the endpoint count and energy equations together.
Net pair losses subtract two dark particles per event. The assigned
reservoir heat is exactly the negative dark energy increment.

For the fixed reservoir temperature T, use L=rho/T−S. Its derivatives
in dark number/energy coordinates are eta and 1/T−1/T_D. At the
endpoint, the gradient pairing of the step equals minus dt times the
direct event entropy production. Convexity of L implies
L(old)−L(new) >= dt*entropy(new); the code checks both ledgers
independently on the numerical steps.

Five initial states and three timesteps per state pass on two positive
quadrature orders. The dark-state susceptibility is also checked
by finite differences. Every displayed step decreases L strictly,
with positive collision entropy and a nonnegative convex remainder.
This is a local fixed-reservoir numerical update, not an expanding
cosmology, evolving visible bath, or arbitrary-step existence claim.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_implicit_reservoir_reference.py
```
