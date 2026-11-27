# Propagate the tiny scalar pole direction without subtracting cold rates

The [attribution check](INDEPENDENT-EARLY-POLE-PAIRED-WARM-FEEDBACK-ATTRIBUTION-2026-11-27.md) found that direct subtraction of cold rates in the short distributed-shell run assigned the wrong sign to the small pole term. First anchor the active gap to the baseline gap by taking the exact difference of the two gap equations, using explicit count increments and rationalizing the difference of inverse one-particle energies. Anchoring agrees with the independently solved masses within 2.2e−19 in these runs but does not remove the spurious rate difference. The dominant loss of precision here is in reevaluation and subtraction of the rate, not just in independently solving the two gaps.

In a second, opt-in experiment, evaluate the cold decay packets at the common baseline pole and add the active pole correction using a centered tangent, multiplied by the gap increment. Apply the tangent to both the cold scalar count source and the reciprocal hard moment deposition, adding it to the differential before combining the full states. The warm packets still use the exact anchored active gap. This is a controlled near-null linearization of the cold term, not a replacement for the ordinary kernel at finite pole shifts.

| Native hard snapshot, six joined steps, t=2e−6 | Anchored direct cold feedback | With reciprocal cold tangent |
| --- | ---: | ---: |
| q=.1–5 markers | −7.40391e−19 | −7.41416e−19 |
| Same shells, distributed q4 | −5.56534e−19 | −7.10356e−19 |

Both tangent histories match the previous independent local-slope estimate to the precision displayed. The magnitude changes by about 4.2% between the representations in this diagnostic, not the roughly 25% implied by raw direct subtraction. Doubling the tangent probe changes the distributed feedback by less than .01%. The terminal comoving energy defects are within 3e−16. This test does not establish a general precision bound or an observational effect.

From `kinetics/`, with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`:
```sh
p=early_pole_live_scalar_exterior_paired_difference_reference.py
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
python -W error "$p" "$a" --joined --duration 2e-6 --steps 6 --order 8 --anchor-pole --linearize-cold-pole
python -W error "$p" "$a" --joined --duration 2e-6 --steps 6 --order 8 --anchor-pole --linearize-cold-pole --warm-order 4
python -W error "$p" "$a" --joined --duration 2e-6 --steps 6 --order 8 --anchor-pole --attribute-feedback --stable-pole-tangent --warm-order 4
python -W error "$p" "$a" --joined --duration 2e-6 --steps 6 --order 8 --anchor-pole --linearize-cold-pole --warm-order 4 --cold-pole-relative-probe 1e-4
```

Raw: `results/early-pole-paired-cold-pole-tangent-2026-11-27.json`.
