# Resolved momentum in the finite cooling-history reference

This is an independent follow-up to the
[finite cooling-clock calculation](INDEPENDENT-FINITE-COOLING-REHEAT-PORTAL-2026-11-26.md).
It keeps the same vacuum lepton annihilation, masses, visible FD factors,
fixed-count normalization and radiation table. The 100 MeV parent decays
to an 86 eV scalar and a stable 20 MeV fermion; direct light-pair
production also adds stable fermions. For the finite branch, temperature
scales as $a^{-3/8}$ while $H\propto a^{-3/2}$. No heating branch,
scattering, inverse, or scalar-pair production is included.

## Count deposition and streaming

At each production node the positive angular density is $N+A P_2(x)$
for a prompt produced stable pair member, with lab energy $C+Dx$.
A boosted two-body decay is also affine in its own independent angular
cosine. Clip each momentum-bin edge analytically in angular space
after redshifting from its production epoch. Three Gauss nodes exactly
integrate both the count and squared momentum on each clipped direct
interval, since the integrands have degrees two and four. The heavy
production boost is integrated by ordinary angular quadrature.
The subsequent decay is uniform in lab energy over each boost interval.
Its two partially covered bins are integrated directly. Full interior
bins receive only nonnegative contributions, accumulated over
half-open ranges before applying the exact energy primitive.
An independent dense angular-clipping test verifies these contributions
at seven temperatures for both species. The main calculation compares
two production/temperature resolutions, with maximum absolute change
of 1.45×10⁻⁴ in the reported CDF probes.

For each positive bin with lower/upper squared momenta $x_0,x_1$,
the geodesic distance from 1 MeV to equality is monotone and concave
in $x=p^2$: the instantaneous speed has the form
$\sqrt{x/(x+C(a))}$, whose second derivative is negative for $x>0$.
Taking the chord between $x_0$ and $x_1$ gives a lower bound, and
the same function at the bin's mean $x$ gives an upper bound. The
overflow lower bound uses its lower edge. Number-weighted sums bound
the represented mean in the prescribed background, tighter than just
taking the global RMS. Quantile brackets use the corresponding bin
edges. These are **distances, not a transfer function**, gravitational
clustering, or likelihood constraints.

Vacuum lifetime is bounded independently in the [prompt-delay audit](INDEPENDENT-FINITE-COOLING-PROMPT-DELAY-2026-11-26.md).

## Results at the 1 MeV snapshot

| Temperature history | Species | Fraction with p/T < 1 | Mean distance to equality |
| --- | --- | ---: | ---: |
| Instantaneous | 86 eV scalar | 5.0×10⁻¹¹ | 19.926–20.050 Mpc |
| Instantaneous | 20 MeV stable fermion | 0.02449 | 0.21906–0.22055 kpc |
| Finite cooling, Tmax/TR=5 | 86 eV scalar | 0.91981 | 2.41460–2.43023 Mpc |
| Finite cooling, Tmax/TR=5 | 20 MeV stable fermion | 0.14356 | 0.15185–0.15289 kpc |

The cold scalar fraction is a feature of dilution and redshifting
during the modeled cooling branch. Over 99.99999% of represented
scalars in the second history are produced before radiation domination.
Their present rest-mass abundance in this omission-controlled model
is only 2.233×10⁻⁹ in Ωh², against 0.120 stable. Thus the colder
scalar spectrum is **not** a dominant scalar dark matter scenario.
These brackets cover the represented streaming mean, not continuum
discretization, omitted production, clustering or observational
uncertainties. Tiny reconstructed low-momentum fractions should
not be treated as relative-precision predictions.

## A positive numerical export, not a transfer calculation

The optional exporter puts each represented bin at its RMS momentum
and widens it with one of two positive lognormal kernels. A median
factor exp(−σ²) preserves the second moment exactly, with σ=0.10 or
0.14, and the table is renormalized to its recorded total number.
Both species and histories close number and p² to better than 10⁻⁷.
Higher moments and within-bin shape are not fixed; e.g. the fourth
moment of these representatives increases by exp(4σ²).
The checked-in manifest records reproducible hashes and tail fractions;
the generated tables are intentionally not checked in.

A preliminary uniform native finite-q grid for the colder scalar
fails the discrete dilation check at 128 nodes (up to about 0.30
in fractional error for the σ=0.10 table). Do not take power from that
grid as a transfer result. A future perturbation calculation needs
independent grid and shape controls and a treatment of the genuinely
tiny abundance at its production value.

## Reproduction

Use the same EOS input and JSON as the previous finite-clock calculation.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_reheat_distribution_reference.py \
    EOS.dat results/finite-cooling-reheat-portal-2026-11-26.json
OPENBLAS_NUM_THREADS=1 python -W error test_finite_cooling_reheat_distribution_reference.py EOS.dat
OPENBLAS_NUM_THREADS=1 python -W error export_finite_cooling_class_psd_reference.py EOS.dat results/finite-cooling-reheat-distribution-2026-11-26.json psd
~~~
