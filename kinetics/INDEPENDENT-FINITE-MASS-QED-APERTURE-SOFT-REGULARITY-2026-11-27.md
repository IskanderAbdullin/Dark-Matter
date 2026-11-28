# A stepped Born aperture is smooth away from tuned occupation edges

The [leading-soft subtraction](INDEPENDENT-FINITE-MASS-QED-SOFT-SUBTRACTION-2026-11-27.md)
requires the thermal symmetric second difference of the measurement
to decay sufficiently rapidly at zero soft shift. A discontinuous
occupation in one shell does not by itself decide the regularity
of the *integrated* aperture.

At fixed scalar momentum Q=4.1e−6, define an explicit test observable:
translate only the scalar energy E_s by a signed amount and
recompute the exact vacuum Born support and the integral of f_P−f_b
over parent energy. The induced invariant scalar mass is updated
so the support inversion is internally consistent. All unchanged
parent and shifted daughter shell boundaries are integrated exactly
as piecewise-constant intervals with 77-digit arithmetic. This is
a diagnostic energy translation, not the actual angle-dependent
recoil map of either radiative process.

* At the two native archive poles, the symmetric difference is
  quadratic. Dividing by the squared shift stabilizes at about
  2.18e5 (fine) and 5.94e5 (third) in the specified test units,
  between shifts 1e−13 and 1e−15.
* Aligning the lower kinematic boundary exactly with an upward
  parent-bin jump changes that second difference to linear in
  absolute shift. Its coefficient follows from the occupation
  jump times the lower endpoint's derivative. The two tuned
  coefficients −59.83 and −104.27 agree with the independently
  integrated aperture to 4.8e−12 relative at shift 1e−15.
  The required energy detunings from the native poles are
  +4.46e−6 and −1.54e−6; this is not an observed native cusp.
* Averaging the tuned aperture over a finite top-hat energy
  width produces finite curvature, given by the difference
  of the aperture slopes at the two averaging boundaries
  divided by the width. Its magnitude scales inversely with
  the width as the bin is narrowed. This is a defined
  measurement regularization, not physical photon damping.

The example distinguishes smooth, stepped, and integrated
measurements before attempting soft matching. It does not certify
infrared safety for the full two crossed recoil maps or for the
nonequilibrium radiative kernel.

Replay with the previously published fine and third archives:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_aperture_soft_regularity_reference.py
```
[Raw receipt](results/finite-mass-qed-aperture-soft-regularity-2026-11-27.json).
