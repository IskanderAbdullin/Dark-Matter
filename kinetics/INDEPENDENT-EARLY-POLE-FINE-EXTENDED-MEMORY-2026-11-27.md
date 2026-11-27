# Five-channel fine history at triple the exposure

Extend the [matched 72-step comparison](INDEPENDENT-EARLY-POLE-JOINED-LOCAL-SPECTRAL-MEMORY-2026-11-27.md) from 2.4e−5 to 7.2e−5 at the same macro width, 3.33333e−7. Both use the same fine hard momentum mesh, fixed scale, prepared inputs, finite channel supports, and unlike response refresh. Each of the first 73 archived count values, including the initial count, is exactly equal to the corresponding value of its previous 72-step branch.

| Quantity at 7.2e−5 | Archived steps | Count-conserving continuous |
| --- | ---: | ---: |
| Integrated scalar gain | 7.845251615e−6 | 7.833774694e−6 |
| Final interval mean rate | .117717978 | .117535156 |
| Final instantaneous rate | .117732426 | .117549529 |

The integrated difference is 0.1463% of the step result, versus 0.1835% at the shorter endpoint with the same width. The final interval-rate gap instead increases from about 0.069% to 0.155%. Final parent and daughter relative occupation L1 differences in the comparison aperture [1.39,1.48] increase from 0.087% and 0.164% to 0.193% and 0.402%. The final full common-support gaps reach 0.731% and 0.546%, respectively; the photon gap is 0.0291%. Maximum absolute energy and entropy ledger residuals remain below 2.5e−13 and 2.0e−12.

Both finite-band rates remain positive and rising near the last sample. Close integrated gains neither imply erased local spectral differences nor show freeze-out. This is a longer calculation at the same step size, not a time refinement or a continuum, in-medium, or expanding calculation.

From `kinetics/`, replay:
```sh
for p in steps cell-conservative; do
  OPENBLAS_NUM_THREADS=1 python early_pole_cold_screened_split_reference.py fine \
    --duration 7.2e-5 --count 216 --same-soft-strength 1 --compton-strength 1 \
    --conversion-strength 1 --screening-update full --pair --projection "$p" --snapshots > "$p.json"
done
OPENBLAS_NUM_THREADS=1 python early_pole_joined_history_memory_reference.py steps.json cell-conservative.json
OPENBLAS_NUM_THREADS=1 python early_pole_joined_snapshot_aperture_reference.py steps.json cell-conservative.json
```

The bundled full histories and both derived diagnostics are in `results/early-pole-joined-fine-216-memory-2026-11-27.json.gz`.
