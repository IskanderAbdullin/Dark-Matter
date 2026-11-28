# Finite-frequency integral of the signed bare cold witness

The [pointwise cold witness](INDEPENDENT-FINITE-MASS-QED-SIGNED-COLD-WITNESS-2026-11-28.md) compares a fixed scalar mode q=4.1e−6 and two parent distributions on the same refined boundaries. Keep those inputs and the bare real leading-soft subtraction, and now integrate over the explicit photon-frequency band [3e−6,1e−5] using dω/ω. The photon population is the same trial Bose factor at Tγ=.9 as in the pointwise test. The .85 mixture is along the same bounded, moment-preserving reconstruction path; it is not an inferred time-dependent parent distribution.

Integrate the spontaneous, stimulated and combined parts independently, and split at each physical angular or radial support knot used by the retained kernel. All panels are recorded at two distinct parent, cosine, azimuth and frequency orders. The mesh changes in the band-integrated combined value are at most 3.19e−10 absolute. The same common-temperature directional-factor control is recorded on each panel; it is not a claim of detailed balance for the archived states.

| Archive | Frequency panels | Native normalized remainder | .85 mixture | Interior zero fraction |
| --- | ---: | ---: | ---: | ---: |
| Fine | 61 | −9.94436e−5 | +8.80613e−3 | .00949148 |
| Third | 125 | +3.02722e−5 | +2.32603e−3 | None |

Every refined native panel is negative on fine and positive on third; the .85 mixture is positive in all panels of both archives. Affinity in the parent occupation on identical refined support fixes the path zero on fine, with a change of 3.00e−8 in mixture fraction between the two quadratures. The native sign difference is therefore not explained by the measured frequency-quadrature error. Neither a global sign nor a crossing independent of the archive or parent path has been established.

This band lies below the plasma gap of the other frozen-pole controls. It is a normalized, signed **bare real** remainder with the leading subtraction only, before the coupling and scalar occupation factors. It supplies neither a propagation law below the gap nor a physical rate there. It does not add the virtual term, match sequential overlap, or fix the sign of an inclusive correction.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error finite_mass_qed_signed_cold_witness_band_receipt_pilot.py --packed
# Regeneration (each panel file is restartable):
for g in fine third; do
  for m in 0 1; do
    python -W error finite_mass_qed_signed_cold_witness_band_pilot.py \
      "$g" "$m" "witness-$g-band-mesh$m.jsonl"
  done
done
python -W error finite_mass_qed_signed_cold_witness_band_receipt_pilot.py .
```
The [receipt](results/finite-mass-qed-signed-cold-witness-band-2026-11-28.json) replays from four exact compressed panel logs.
