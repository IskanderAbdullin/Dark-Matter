# Check the cold proxy for the archived stable fermion

The [matched native scalar comparison](INDEPENDENT-WARM-DIRECT-NATIVE-EXACT-K-2026-11-26.md)
uses CDM for the dominant stable species. Test that proxy
separately with the represented 20 MeV stable-fermion phase
of the corrected 7.394 MeV archive. This does not generate
another full Boltzmann metric or write a new perturbation
history for the original portal model.

## Homogeneous replacement

At equal present density subtract the baseline CDM energy
using the positive form of sqrt(1+x²)−1. At the starting
archive the stable phase makes H larger by 2.5802e−8 in
the prescribed radiation background. At equality this
difference is only 1.2663e−15, and the accumulated conformal
time difference from the archive to equality is −5.59e−12
Mpc. Two original-cell quadratures reproduce the results.
No velocities before the archive or unresolved phase are
included in this replacement.

## Forced tracer against its own cold limit

Use the same recorded metric for the finite-mass stable
tracer and the analytic q/m→0 limit of its collisionless
operator. Start with unit energy contrast and zero peculiar
velocity for both. The stable initialization divides its
dilational derivative by 3+<v²>_energy so the initial energy
perturbations are exactly matched. The cold limit keeps
the finite product of j1(k chi) and 1/v in the source
rather than setting them independently to zero. This leaves
a source containing a times the conformal integral of
1/a. Compare at z=20; use a second forcing-only control
with zero initial contrast.

| Prescribed fixture k (Mpc⁻¹) | Fractional stable-minus-cold tracer amplitude |
| --- | ---: |
| 1.3472 (2 h/Mpc) | −2.62e−11 |
| 2.6944 | −1.30e−10 |
| 6.736 | −1.02e−9 |

The rms ballistic path to equality after the archive is
.2369 kpc. The independent initial-only spherical bound
1−j0(k chi)≤k² chi²/6 gives ≤1.70e−8 at 1.3472 Mpc⁻¹;
it is not a bound on the forced tracer. The matching
forced comparison is smaller. These explicit phase,
initial-data and background tests support the chosen CDM
proxy in this setting, not a certified correction to
the native total matter power or a bound for different
stable production spectra.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
A=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
python -W error warm_native_stable_background_proxy_reference.py EOS.dat "$A"
python -W error warm_native_stable_cdm_limit_reference.py EOS.dat "$A" results/los-metric-off/*.npz
# The CDM-only fixture is generated in the recorded metric survey.
python -W error warm_native_stable_cdm_limit_reference.py EOS.dat "$A" --cdm-driver \
 cdm_native/cdm_los_00_perturbations_k1_s.dat
~~~

[Homogeneous replacement](results/warm-native-stable-background-proxy-2026-11-26.json);
[two recorded fixtures](results/warm-native-stable-cold-driver-2026-11-26.json);
[separate CDM-only fixture](results/warm-native-stable-cold-k2-driver-2026-11-26.json).
