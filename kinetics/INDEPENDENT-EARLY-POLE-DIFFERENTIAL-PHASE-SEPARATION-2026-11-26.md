# Equal charges cancel a common phase, but the paths can separate

At the resonant parent energy midpoint, on the same pole dispersion as the
finite-aperture check, the exact velocity vectors differ by 6.0361e−4 at
illustrative coupling g=.005. The angular detuning half-span is 4.10e−6 T.
Straight paths separate by 147.22 / T in one inverse half-span. Taking the
inverse ideal massless two-Dirac Debye scale only as a comparison length,
this separation is .601 times that length (.120 at g=.001 and 1.202 at .01).
These are kinematic values, not a calculated decorrelation time.

An explicit imposed Gaussian process with time-white scalar potential and
spatial covariance C(r/L) makes the trajectory effect testable. If the
single phase averages to exp(−Γt), the neutral relative phase averages to
exp[−2Γ times the time integral of (1−C(δv t/L))], by the covariance of
the phase difference. Its decay vanishes on identical paths at any Γ.
A smooth spatial Gaussian gives an early cubic time dependence of the
log-amplitude; a cusp exponential C produces quadratic dependence instead.
Neither case is a stationary spectral ansatz: both paths start together and
their separation changes with elapsed time.

At g=.005 and one inverse angular half-span, the accumulated pair exponent
is .0571 (Gaussian) or .2484 (exponential) of two independent line exponents.
For Γ equal to the angular half-span, the first pair e-fold takes 4.39 or
2.96 times the e-fold time of independent lines. The exact covariance
integrals agree with independent scalar quadrature to 3e−16 absolute.
The two toy choices demonstrate why a single-line exponent and one
screening length alone do not fix the neutral phase response.

No Gaussian-process parameters are fitted to the actual gauge plasma.
The ideal Debye label does not describe magnetic, dynamic, nonequilibrium
or collisionally bent paths. Real gauge dynamics must account for the
current difference, thermal fluctuation spectrum and corresponding vertex
contributions. No cosmological production rate is derived here.

Reproduce from kinetics/ with early_pole_differential_phase_separation_reference.py.
Raw output: results/early-pole-differential-phase-separation-2026-11-26.json.gz.
