# Late rate agreement and residual charged-spectrum differences

In the [existing finite cold window](INDEPENDENT-EARLY-POLE-JOINED-INPUT-MEMORY-2026-11-26.md), the fine-archive comparison is extended to reference exposure 2.4e−5, using 72 equal macro steps. The same two reconstructions are also tested independently on the third archive using 24 macro steps up to 8e−6. Both are fixed-scale runs with photon pair reactions, hard Compton and hard flavor conversion, unlike screened elastic, and the same-flavor screened pole. Only the unlike soft polarization is actively refreshed. The input keeps the same scalar counts in the two projections of any given archive, but the two native archives have different initial states and grids.

| Archive, exposure | Integrated scalar gain: step / continuous | Final interval rate: step / continuous |
| --- | --- | --- |
| Fine, 2.4e−5 | 2.3574005e−6 / 2.3530754e−6 | .10873799 / .10866306 |
| Third, 8e−6 | 7.1152660e−7 / 7.1413523e−7 | .09716174 / .09718144 |

The final interval's relative rate difference is 0.069% on fine and 0.020% on third. The fine interval rates stay within 0.1% after exposure 1e−6; both final rates are higher than their respective initial rates. On fine at the final snapshot, the occupation relative L1 in the common momentum aperture [1.39,1.48] is 0.087% for the parent and 0.164% for the daughter, down from 6.06% and 19.82% at initialization. Over all common live support, however, the final relative L1 is still 0.680% and 0.410%. Initially local occupations relax rapidly in these two trajectories, but the observable's agreement does not imply identical hard spectra or establish generic relaxation of all inputs. In particular, a [separate frozen bare witness](INDEPENDENT-EARLY-POLE-NATIVE-APERTURE-WITNESS-2026-11-27.md) demonstrates that matching three archived moments does not force the initial source sign.

The independent live-aperture replay uses the variational dressed energy gradients and scalar Hartree gap already used in the joined dynamics. It reconstructs the separate forward occupation, inverse occupation and spontaneous integrals instead of inferring spectral agreement from a scalar difference. Per-mode instantaneous sources agree with the saved derivatives within 1.4e−17 absolute. For the fine step branch, the separate stimulated contributions (parent, daughter magnitude) move from (.112542, .033906) to (.131802, .022998). The final step-minus-continuous differences in these two contributions are +.00011453 and +.00003865, whose difference explains the +.00007589 final instantaneous net source difference. This is a decomposition of the modeled collision operator, not an independent physical rate calibration.

For third at exposure 8e−6, corresponding no-elastic control gains are 5.7714265e−7 and 6.8591342e−7; adding both elastic sectors shifts them by +23.28% and +4.11%, respectively. A separate fine native-step no-elastic control at 2e−6 yields 1.5710013e−7 against 1.5212960e−7 with elastic (−3.16%); at 8e−6 the published matched comparison gives +10.72%. The sign is not constant across these fixed-clock finite-observation windows, and this comparison does not identify the crossing time. We do not infer temporal convergence at 2.4e−5 from a 72-step result alone. All quoted tested sources remain positive at sampled late times and there is no freeze-out inferred here. The maximum measured energy/entropy residuals of the five-channel long fine and third histories are below 2.5e−13 / 2.2e−12.

From `kinetics/`, with one BLAS thread, repeat for `cell-conservative` and change `GRID T N` independently between `fine 2.4e-5 72` and `third 8e-6 24`:

```sh
GRID=fine; T=2.4e-5; N=72; P=steps
OPENBLAS_NUM_THREADS=1 python early_pole_cold_screened_split_reference.py "$GRID" \
  --duration "$T" --count "$N" --same-soft-strength 1 --compton-strength 1 \
  --conversion-strength 1 --screening-update full --pair --projection "$P" --snapshots > "$P.json"
OPENBLAS_NUM_THREADS=1 python early_pole_joined_history_memory_reference.py steps.json cell-conservative.json
OPENBLAS_NUM_THREADS=1 python early_pole_joined_snapshot_aperture_reference.py steps.json cell-conservative.json
# Matched 8e-6 third no-elastic control on both projections:
OPENBLAS_NUM_THREADS=1 python early_pole_cold_screened_split_reference.py third \
  --duration 8e-6 --count 4 --strength 0 --compton-strength 1 --conversion-strength 1 \
  --pair --projection "$P"
# Fine 2e-6 no-elastic control: same command with the following three substitutions.
# third -> fine; --duration 8e-6 -> --duration 2e-6; --count 4 -> --count 2. Use steps.
```

Raw: `results/early-pole-joined-fine-72-memory-2026-11-27.json.gz` and `results/early-pole-joined-third-24-memory-2026-11-27.json.gz`. The bundled 24-step third history predates collecting optional occupations, so its comparator uses scalar history alone.
