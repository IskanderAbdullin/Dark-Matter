# Separate direct stable and scalar effects on a fixed metric

The [shifted scalar result](INDEPENDENT-FINITE-COOLING-SHIFTED-SPECTRUM-2026-11-26.md)
is a density contrast for a trace population. It need not be the
largest direct kinematic change in the same composition. Here propagate
the 20 MeV stable distribution as an additional collisionless tracer on
the identical CDM-only metric. Use the shifted prompt stable shape
for this diagnostic, independent of its overall normalization. At
initial (k\tau\) about 0.07 and 0.15, separately match the native CDM
density monopole and bulk velocity, with no higher initial moments.
For the latter use the linear local boost of an isotropic occupation
((\epsilon/q\) times the momentum derivative). This is not a solution
for stable perturbations before the first match.

## Resolve the small paired difference

Integrate a second, exactly pressureless tracer driven by the **same**
metric, with `delta' = -theta + 3 phi'` and
`(a theta)' = a k² psi`. The pressureless result reproduces
the native CDM density to about 1.3×10⁻⁷ at the tested endpoints.
Resolve the much smaller kinetic difference using the zero-velocity
limit as a control variate inside the kinetic integral:
`j0(x)-1`, `j1(x)/x-1/3`, and the difference of inverse energies
are computed without subtracting nearly equal rounded quantities.
For the cold part the weighted derivative identity is integrated
analytically. Ignoring the shared numerical trajectory and directly
subtracting two independently sampled CDM-size values does not
reliably resolve such a small response.

Both momentum/time refinements and both matched starts were tested.
At z=20 and the three sampled modes, the relative stable suppression
in density contrast is between 1.3×10⁻¹⁰ and 1.02×10⁻⁸.
At the largest tested mode, doubling the rest mass at the same
frozen momentum shape nearly divides the paired correction by four
(the nonrelativistic velocity-squared scaling control). This mass
doubling is a numerical test, not a second production scenario.

## Direct test-particle composition, not a transfer calculation

Hold the dark and baryon density weights of the driver fixed.
Combine the stable tracer and scalar tracer at their own tiny
scalar fraction within the stipulated dark component, and reuse
the baryon trace. The scalar fractions are 2.54×10⁻⁹
(instantaneous) and 2.29×10⁻⁸ (finite). These are rest-count
fractions, not fractions of an inferred total gravitating mass.
At the earlier matched start, the direct squared-contrast change
relative to the pressureless-dark comparator is:

| k [h/Mpc] | Instantaneous | Finite, Tmax/TR=5 |
| ---: | ---: | ---: |
| 4 | −3.20×10⁻⁹ | −4.89×10⁻⁹ |
| 10 | −6.99×10⁻⁹ | −1.12×10⁻⁸ |
| 22 | −2.15×10⁻⁸ | −2.49×10⁻⁸ |

The alternate initial match preserves this scale of change.
These values are **not** observational power predictions, an
independent combined CLASS solution, a bound at unsampled modes,
or a limit on a more abundant or differently produced scalar.
Neither test tracer changes the metric or feeds back onto baryons.

## Reproduction

Use the archived shifted prompt reconstruction at width 0.14 and
the previous CDM-only native driver. The two stable tables are
in that export's prompt directory:

~~~sh
export OPENBLAS_NUM_THREADS=1
python -W error finite_cooling_stable_paired_los_reference.py EOS.dat \
  prompt/{instant,ratio5}-stable-sigma0.14.dat cdm_native/cdm_los_00_perturbations_k{2,4,6}_s.dat >stable.json
python -W error finite_cooling_stable_paired_los_reference.py EOS.dat \
  prompt/{instant,ratio5}-stable-sigma0.14.dat cdm_native/cdm_los_00_perturbations_k6_s.dat --mass-factor 2 >double.json
D=results
python -W error finite_cooling_direct_matter_mix_reference.py stable.json double.json \
  "$D"/finite-cooling-shifted-spectrum-{prompt,mix}-2026-11-26.json \
  "$D"/finite-cooling-shifted-count-audit-2026-11-26.json "$D"/finite-cooling-cdm-los-input-2026-11-26.json \
  "$D"/finite-cooling-shifted-spectrum-{instant,ratio5}-2026-11-26.json
~~~

## Separate high-mode control

Use k=34 and 40 h/Mpc from the separate CDM-only driver
prepared with `--modes 26 30 34 40`. Repeat both projections
on those exact files; do not interpolate between the survey
and high-driver modes. At the earlier start, direct squared-
contrast shifts at k=34,40 are respectively −4.78×10⁻⁸ and
−6.58×10⁻⁸ (instantaneous) versus −4.32×10⁻⁸ and −5.46×10⁻⁸
(finite). The alternate start preserves this ordering. Unlike
the first three modes above, the instantaneous case now
gives the larger direct suppression of the two, because its
20 MeV tracer is faster. This is an ordering of sampled
points, not a continuously located crossover or a prediction
for total matter in a self-gravitating model.

For the high comparison select the new high native files
`cdm_los_00_perturbations_k2_s.dat` and `_k3_s.dat` in both
projections. Use the matching high-driver configuration
`finite-cooling-cdm-los-high-input-2026-11-26.json` and the
two new scalar records in the composition command. Keep
`double.json` as the k=22 mass control and add
`--control-20 stable.json` pointing to the corresponding
k=22 value from the first run.
