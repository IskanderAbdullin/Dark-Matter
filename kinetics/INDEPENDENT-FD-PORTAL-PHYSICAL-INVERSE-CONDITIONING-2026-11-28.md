# Conditioning of the actual local FD number/heat inverse

At each sampled visible temperature and force, compute the true bath-dependent
rate Jacobian J. With s_i=sqrt(J_ii), use the normalized tangent
diag(1/s) J diag(1/s). Its singular values refer to this particular
force/rate scaling, not to a temperature-independent norm.

Prescribe equal-scale small local force perturbations along the weak and
strong right-singular directions. Evaluate the true finite rate targets
on both sides, invert them independently, and compare their antisymmetric
rate differences with the singular-mode prediction. The two positive
quadratures independently reproduce the normalized condition number.

| Fixed visible T (GeV) | Normalized condition number |
|---:|---:|
| .05 | 43.19 |
| .30 | 19.12 |
| .64135 | 492.48 |
| .6882487043 | 419.28 |

In this numerical metric the higher-temperature samples distinguish one
certain force combination poorly relative to the strong direction, even
though the actual local inverse succeeds from a nearby starting state.
These are local sensitivity controls in the finite vacuum lepton model,
not observational uncertainties or bounds on an evolving cosmology.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_physical_inverse_conditioning_reference.py
```
