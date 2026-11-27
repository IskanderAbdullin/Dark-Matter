# Match both early radiation and late kinetic momentum, then test native power

The [conditional native pilot](INDEPENDENT-WARM-DIRECT-NATIVE-PILOT-2026-11-26.md)
showed that preserving just q² does not preserve the
scalar power response. The early relativistic energy is
controlled by mean q, while the leading late kinetic
energy is controlled by mean q². Construct a positive
mixture of two thermal-shaped templates that matches
the direct distribution in both. For equal particle
counts in the two templates, choose their mean
temperature so the first moment matches and their
temperature variance so the second moment matches.
This gives temperatures .33839 and .68310 in photon
units, with mean q=1.60954 and q²=3.75972. These are
numerical scalar spectrum templates, not assumptions
about its quantum-statistical microphysics.

The single q²-matched template has mean q=1.69872.
At a=5e−6, close to the transition, the scalar pressure
ratios for direct, single and two-temperature shapes
are .149267, .156357 and .148997. The recovered
adiabatic sound speed is independently checked by
finite differencing pressure and energy in log a
(agreement better than 1.6e−10 relative). At late times
matching q² controls the leading pressure but leaves
different q⁴ and other distribution information.

| Method at amplified f=1e−3, z=5 | Single temperature minus direct, Delta ln P at k=1 | Two temperatures minus direct |
| --- | ---: | ---: |
| Standard NDF | −2.313e−4 | −3.190e−5 |
| Early-start RK, exact k | −2.124e−4 | −1.506e−5 |

The second line uses the [common early-start control](INDEPENDENT-WARM-DIRECT-NATIVE-EXACT-K-2026-11-26.md),
so it avoids the known large numerical offset from
different starting criteria. The two-moment template
reduces the residual by about 93% on that particular
mode but does not eliminate it. The small-mode
differences are sensitive to precision/choice of
physical and numerical comparisons, so do not infer a
sign at the lowest modes. These are deliberately
amplified abundance tests, not a linear extrapolation
or a prediction of the actual scalar abundance.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_direct_vacuum_two_moment_native_control_reference.py EOS.dat native
/path/to/class native/std-p12.ini >native/std-p12.log
python -W error warm_direct_vacuum_native_early_start_reference.py focus early --opacity 2e-7 --templates native
for x in direct3 thermal3 p12; do
 /path/to/class "early/$x.ini" >"early/$x.log"
done
python -W error warm_direct_vacuum_native_moment_match_audit_reference.py native early
~~~

[Template moments and pressure](results/warm-direct-vacuum-native-two-moment-template-2026-11-26.json);
[native comparison](results/warm-direct-vacuum-native-two-moment-power-2026-11-26.json);
[exact-mode native tables](results/warm-direct-vacuum-native-two-moment-raw-2026-11-26.json.gz).

[Unamplified equal-fraction shape check](INDEPENDENT-WARM-DIRECT-NATIVE-PHYSICAL-SHAPE-2026-11-26.md)
uses one matched trace for each of the two spectrum tables.
