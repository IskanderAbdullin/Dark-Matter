# Verify the physical scalar response at exact native wavenumbers

The [first native pilot](INDEPENDENT-WARM-DIRECT-NATIVE-PILOT-2026-11-26.md)
used NDF and an interpolated common grid. In an independent
narrower run, explicitly request k=.002,.01,.03,.1,.3,1
h/Mpc. All six appear exactly in every native power file.
Use explicit RK, 96 momenta, lmax=24, tolerance 9e−10 and
the identical conditional direct warm scalar PSD and
abundance. No k interpolation is needed for this audit.

The trace-species offset at k=.002 remains about
−1.16e−5 against the model without the species. This
rules out the common-grid interpolation as its sole cause.
The recorded first perturbation output at this mode is
tau≈56.87 Mpc without the species, versus .087697 Mpc
for trace and physical fractions. The trace and physical
models start at essentially the same scale factor
1.89213e−7 on every saved mode. CLASS has a separate
relativistic-initialization requirement for an NCDM
component, even one carrying negligible density; the
observed starts reinforce using the matched trace rather
than an unmatched zero-species baseline.

| k (h/Mpc) | Exact physical minus trace Delta ln P, z=5 | Same at z=20 |
| --- | ---: | ---: |
| .002 | −2.48e−9 | −2.49e−9 |
| .1 | −7.97e−8 | −7.94e−8 |
| .3 | −9.52e−7 | −9.49e−7 |
| 1 | −5.088e−6 | −5.064e−6 |

Compared to the earlier standard pilot, the largest
absolute difference is 3.61e−8; at k=1 this is only
about 0.71% of the signal. Centering the amplified
1e−4 case on the exact-mode trace predicts both
sets within 7.31e−9 absolute. The sign at k=.002
flips relative to the standard pilot and is not resolved.
The checks support a few-parts-per-million decrease at
k=1 for this prescribed scalar spectrum, not a revised
production or observed matter constraint. The narrowed
run stops below k=2 and does not certify the higher
wavenumbers or the separate larger hierarchy.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p focus
python -W error warm_direct_vacuum_native_focus_reference.py native focus > focus/spec.json
for name in base trace physical amp4; do
 /path/to/class "focus/$name.ini" >"focus/$name.log"
done
python -W error warm_direct_vacuum_native_focus_audit_reference.py focus \
 results/warm-direct-vacuum-native-class-pilot-2026-11-26.json
~~~

[Exact-mode and starting-time audit](results/warm-direct-vacuum-native-exact-k-focus-2026-11-26.json);
[base and trace tables](results/warm-direct-vacuum-native-focus-raw1-2026-11-26.json.gz);
[physical and amplified tables](results/warm-direct-vacuum-native-focus-raw2-2026-11-26.json.gz).
