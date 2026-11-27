# A matched native pilot resolves the trace-species baseline issue

The [matched setup](INDEPENDENT-WARM-DIRECT-NATIVE-BASELINE-2026-11-26.md)
puts the scalar at its conditional physical fraction
2.099185e−6 of present dark density and replaces the same
amount of CDM. A first native CLASS 3.3.4 calculation uses
synchronous gauge, full non-cold hierarchy up to l=20,
70 momentum nodes, and perturbation tolerance 2e−8.
All quoted modes are in h/Mpc. The reported non-cold
density is 2.51902e−7 in omega units, as requested, and
the native relativistic Delta N_eff=1.9681e−7 lies inside
the independently estimated shell range.

Do not subtract only the zero-species result. Even inserting
a trace fraction 1e−10 changes the solver path: its difference
from the zero-species case reaches 1.16e−5 at k=.002. That
is far larger than the trace fraction and should not be
interpreted as a physical small-fraction response. Use the
trace run with the same phase-space layout, hierarchy and
tolerances as the physical baseline. This experiment isolates
a solver setup/baseline sensitivity; it does not diagnose the
complete CLASS initialization mechanism.

| k (h/Mpc) | Physical minus trace, Delta ln P at z=5 | Same at z=20 |
| --- | ---: | ---: |
| .002 | +1.04e−8 | +1.05e−8 |
| .1 | −5.04e−8 | −5.02e−8 |
| .3 | −9.16e−7 | −9.13e−7 |
| 1 | −5.05e−6 | −5.03e−6 |
| 2 | −7.86e−6 | −7.78e−6 |

The 1e−4 amplified case, centered on the same trace
reference, predicts the physical values with maximum
absolute deviation 5.71e−9 on the sampled modes.
This helps distinguish the higher-k few-parts-per-million
suppression from the solver-presence offset. At the
smallest modes do not interpret the sign as a resolved
physical prediction. A moment-matched FD comparator
at the same amplified 1e−3 fraction gives a different
power response, e.g. lower Delta ln P by about 2.31e−4
at k=1 and z=5 than the direct shell spectrum. Same
squared momentum and rest fraction do not determine
native transfer at those test settings.

These are first-tier linear native results, not completed
momentum, hierarchy, or small-k convergence. In this pilot
the scalar spectrum is already present and collisionless,
no other dark distribution is regenerated, and the CDM
replacement is a conditional proxy. The tighter runs are
kept separate and should replicate (or correct) the result.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# Generate "native" according to the matched setup note; place the
# provided JSON run specification at native/spec.json.
for name in base trace physical amp4 amp3 amp2 thermal3; do
  /path/to/class "native/std-$name.ini" >"native/std-$name.log"
done
python -W error warm_direct_vacuum_native_class_audit_reference.py native /path/to/class --partial
~~~

[Provisional signed native comparison (temporary paths elided in logs)](results/warm-direct-vacuum-native-class-pilot-2026-11-26.json);
[raw base, trace, physical](results/warm-direct-vacuum-native-pilot-raw1-2026-11-26.json.gz);
[raw amplified and thermal outputs](results/warm-direct-vacuum-native-pilot-raw2-2026-11-26.json.gz).

## Standard-tier PSD smoothing control

At the same amplified scalar abundance 1e−3, replacing
the direct spectrum with the held q64 spectrum lowers
Delta ln P by about 4.55e−4 at k=1, z=5. Altering the
positive PSD smoothing width from .07 to .11 in log q
while preserving the chosen squared momentum changes
either individual power response by at most 1.22e−5,
and the held-versus-direct contrast by at most 3.27e−6
over the sample. These are controls on the specified
smoothing choices, not a bound on all reconstructions
or the higher native tier.

[Smoothing comparison](results/warm-direct-vacuum-native-pilot-smoothing-2026-11-26.json)
and [native tables and elided-path logs](results/warm-direct-vacuum-native-pilot-shape-raw-2026-11-26.json.gz).

[Exact native modes and earlier initialization](INDEPENDENT-WARM-DIRECT-NATIVE-EXACT-K-2026-11-26.md)
repeat the key pilot check with RK and no k interpolation.
