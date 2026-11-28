# Homogeneous exponential damping does not improve these forced expansion tests

Compare the unfiltered and filtered closed FD kicks within the same entropy/pressure-work expansion split, the same physical timelike FD ODE, and the same initial temperatures and perturbations. Check the reference at two independent integration tolerances and on two collision meshes. The largest difference of any compared split endpoint between the collision meshes is 3.89e-13 in its reported coordinates. Each method passes its own discrete entropy ledger and comoving energy/work ledger. The filtered entropy production remains the modified production, not the physical finite-step rate.

| Initial visible temperature (GeV) | Eight-step unfiltered scaled state error | Eight-step filtered scaled state error |
| ---: | ---: | ---: |
| 0.64135 | 1.61e-7 | 3.24e-5 |
| 0.68825 | 1.62e-7 | 5.23e-5 |

Thus the filter, which damps the homogeneous fast mode well in the closed [unforced test](FD_PORTAL_TABULATED_EXPONENTIAL_FILTER.md), makes these two moving cases less accurate. The sampled eight-step ratios are approximately 201 and 323. This rules out choosing the filter from the homogeneous test alone; it does not show that the unfiltered method is uniformly preferable under arbitrary stiffness or forcing.

An exact scalar comparison gives a mechanism. For constant forcing b and decay lambda, an additive half-drift, homogeneous collision and additive half-drift maps `u -> R*u + (dt/2)*(1+R)*b`. With `R=exp(-z)` and `z=lambda*dt`, the numerical steady response is `(z/2)*coth(z/2)` times the true `b/lambda`. The unfiltered midpoint Cayley factor preserves the correct steady response in this affine model, though it can change the sign of a stiff transient. At the two eight-step probe values, the scalar factors are 1.249 and 1.806, and are not estimates of the nonlinear errors.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error fd_portal_stiff_expansion_filter_probe.py eos-range.dat
python -W error fd_portal_affine_split_forcing_reference.py results/fd-portal-stiff-expansion-filter-probe-2026-11-28.json
```

This controlled timelike-only study has a short expanding interval and no abundance inference.
