# Why the affine forcing correction cannot be added to an isolated collision kick

The exact compensation to the [constant-forced scalar split](FD_PORTAL_STIFF_EXPANSION_FILTER_PROBE.md) is a constant offset in its full-step map. Write `R=exp(-z)`, `u*=b/lambda` and `z=lambda*dt`. The required offset is `c=((1-R)-z*(1+R)/2)*u*`. It restores the exact affine evolution for every initial amplitude if inserted between the two additive half-drifts.

However, when the isolated middle collision enters at equilibrium `u=0`, the offset sends it to `u=c`. The passive quadratic free energy then increases by `c^2/2>0` for nonzero forcing and positive z. Thus this naïve assignment to the collision substep cannot preserve the original passive entropy property. It is not a recommended physical correction; a forcing-aware method would have to handle expansion and collisions together or otherwise repartition the forcing and audit the conserved budgets.

The executable checks the exact full-step identities over eight z values and both signs of forcing, and reports the obstruction at the z values used in the preceding nonlinear probe. These are mathematical control cases, not estimates of physical entropy error.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_affine_forcing_correction_guard.py \
 results/fd-portal-affine-split-forcing-2026-11-28.json
```
