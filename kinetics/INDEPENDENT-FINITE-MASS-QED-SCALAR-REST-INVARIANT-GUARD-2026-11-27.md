# Replay a finite charged invariant guard in two Lorentz frames

Fix a manufactured internal occupation window with constant weight
between two nearby energies in the neutral scalar rest frame.
Keep the exact same window, scalar integrand and finite imaginary
charged invariant in two independent coordinate descriptions.
This is a geometry test, not the archive occupation or a damping law.

In the scalar rest frame the uncut charged denominator is linear
in the internal radial energy displacement, with slope ±2 m_scalar.
The photon angle can be integrated directly because the chosen
occupation window is isotropic in this frame. In the bath frame
the same occupation becomes a polar strip. Its bounds depend on
the internal radial momentum, and cannot be replaced by an isotropic
bath-frame occupation. Use the exact complex azimuth/polar primitive
on that moving strip, with the same finite charged invariant.

To resolve the lab integral, split at the boosted radial aperture
edges and at the separate points where the narrow photon angular
peak crosses a strip boundary. Dropping the latter breaks leaves
a measurable quadrature error even with many nodes. After resolving
them, independent coordinate results agree in both real and imaginary
parts to 5.2e−11 absolute, and the two lab meshes to 2.4e−11.
The tests cover three finite invariant guards from 1e−6 to 1e−10
and both charged topologies on two fixed events.

Equality at fixed invariant guard is expected from using the same
integral. This does not compare different zero-width paths or allow
a Lorentz-frame-dependent hard cutoff to be assumed equivalent.
No physical spectral width or inclusive thermal rate is inferred.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_scalar_rest_invariant_guard_reference.py
```
[Raw data](results/finite-mass-qed-scalar-rest-invariant-guard-2026-11-27.json).
