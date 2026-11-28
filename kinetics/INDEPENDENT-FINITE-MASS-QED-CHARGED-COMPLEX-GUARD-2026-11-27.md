# Convert the smooth charged-pole guard to a sharp radial guard

The [bare angular/radial reference](INDEPENDENT-FINITE-MASS-QED-CHARGED-CLOSED-RADIAL-PV-2026-11-27.md)
uses a principal value in the uncut charged denominator and a sharp
radial exclusion. Here give the charged denominator a positive
imaginary invariant i epsilon, while keeping the photon bare. This
is a mathematical regulator comparison, not a physical damping width.

After exact azimuthal averaging, write the photon discriminant
as R²[(u−c_hat)²+a²]. Set u=c_hat+a sinh y and v=exp(y).
With z the complex charged root, h=(z−c_hat)/a, the two transformed
roots are h±sqrt(1+h²) and have product −1. The polar pole integral
is a difference of complex logarithms along a positive real v
interval. The imaginary part approaches the polar delta term with
the sign dictated by the positive i epsilon. Direct finite-complex
adaptive quadratures agree to 2.1e−12 relative; real PV and imaginary
delta limits agree to 2.5e−11. Regulating the original charged
denominator or the partial-fraction pole differs by a finite epsilon
term that tends to zero in the tested window.

At a tuned native step, the finite-width radial results reproduce the
leading coefficients of log(1/epsilon): +.72270134 and +.61595816
versus +.72270204 and +.61595881 for the sharp radial prescription.
At the two untuned central events this logarithmic slope vanishes.

To compare finite parts as well, use the local width a≈A|delta| and
real pole shift B delta about the external energy. Define
delta0=epsilon/(2Qp A), and let S(t) be the real local infinite-angle
Hilbert integral at shift B/A−i/t divided by its undamped value.
The conversion constant is

    K = integral_0^infty [S(t)−1_(t>1)] dt/t.

Subtract J[log(L/delta0)+K] from the smooth radial answer when
the local step coefficient is J and radial half-width is L.
Equivalently the effective sharp exclusion is delta0 exp(−K).
On the fine daughter and third parent tuned events K is about 12.048
and 8.874. This finite conversion is necessary even when both raw
slopes match. After conversion, the two finite parts agree within
1.7e−6 absolute; at the two untuned controls no log conversion is
required and the finite limits match at the printed precision.

The result compares two bare prescriptions near fixed archive
events. It does not pick a thermal real-time contour, produce a
dressed spectral function, or turn the archive step into a
physical divergence.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_charged_complex_guard_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-charged-complex-guard-2026-11-27.json).
