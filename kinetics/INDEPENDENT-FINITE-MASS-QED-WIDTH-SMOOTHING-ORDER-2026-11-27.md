# How a continuous transition removes the local real regulator offset

The [correlated two-width experiment](INDEPENDENT-FINITE-MASS-QED-CORRELATED-OVERLAP-GUARDS-2026-11-27.md)
found a finite bare real offset for an exact occupation step. This
test changes only the local occupation profile, not either guard
or the causal signs. It is a model of the crossover and is distinct
from the number-and-energy-conserving archive reconstruction.

In the analytic local scaled model, the real difference of two-width
and charged-only kernels is odd in radial displacement. For a
continuous transition centered on the old step, replace the step
jump in the paired integrand by its value times G(x/r), with
r = smoothing width divided by the charged overlap radial scale.
Compare G=erf(x/r) with G=t²(3−2t), t=min(1,x/r); their central
slopes differ. Normalize by the independently integrated exact
step offset at the same regulator ratio lambda=1.

At very small r the integral recovers the full step offset. At
r=1 it still retains about .96–.98 of that value, depending on
the fixture and shape. At r=10,000 the fraction is about .100–.114
on the fine daughter versus .010–.022 on the third parent.
At r=10^8 all tested fractions are below 3.2e−5. The two
meshes and scaled supports agree within 2e−6 absolute.
This behavior illustrates that keeping the transition width
fixed as the charged regulator vanishes removes the local
real step offset, whereas shrinking the transition with the
regulator can retain a finite, profile-dependent value.

This is a local overlap test with manufactured profiles, not a
physical collision rate or a prescription for a hot-medium width.
It does not remove or bound the separate reconstruction sensitivity
away from the coincidence.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_width_smoothing_order_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-width-smoothing-order-2026-11-27.json).
