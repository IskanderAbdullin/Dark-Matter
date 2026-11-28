# Fine charge refinement: the cold-band magnitude and neutral sign interval move

Compare the independently evolved fine d1 and fine d2 charge charts at x=.05. Leave the scalar shell edges unchanged, but use each chart's own evolved occupations, fitted scalar mass and transverse pole at fixed k=.0015. Both collision integrals retain q in [1e−6,1.5e−5]. The two d2 meshes agree on all nine integrated coefficients within 5.78e−12 relative across 49/49 panels. The matched d1 result has 27 panels and an accepted comparison within 2.18e−11. The saved endpoint states independently reproduce both pole records and exact panel edges.

| Held numerical state | Forward number relative to original fine archive | Conditional number-neutral occupation | Relative energy-minus-number neutral split |
| --- | ---: | ---: | ---: |
| Fine d1 | 6.45841084 | 3.34922719e−9 | +4.47663e−8 |
| Fine d2 | 1.30271123 | 2.90893020e−9 | +5.34301e−8 |

The same-band gross forward coefficient falls **79.829%** from d1 to d2; the inverse coefficient falls **82.481%**. The net number expression is `n_gamma D - B`. For any held occupation between the two number roots, it is negative on d1 and positive on d2. At their constructed geometric midpoint (~3.12132473e−9) the net is −6.80% and +7.30% of the corresponding gross inverse term. The midpoint is a diagnostic selected to illustrate this interval, not a predicted photon occupation.

We also bound where the collision integrals are supported by using only positive panel weights, with no interpolated subpanel shape. The median of gross forward number is enclosed at q=[3.5685,4.7783]e−6 for d1 and [3.3060,4.3026]e−6 for d2; these intervals overlap, so the data do not order the medians. The fraction at q≥1e−5 is exactly a panel sum in both charts: 3.386e−6 versus 4.876e−6 of forward number and 3.514e−5 versus 5.045e−5 of forward excess. Fractional and absolute contributions are different: the total forward coefficient fell substantially. No extrapolation beyond this interval or claim of a spatially converged plasma source follows.

The [d2 receipt](results/early-pole-fine-d2-endpoint-continuum-cold-band-2026-11-28.json) and [matched contrast](results/early-pole-fine-d2-cold-panel-contrast-2026-11-28.json) retain exact values. From `kinetics/` with the endpoint dependencies installed:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_cold_band_panel_contrast
gzip -dc results/early-pole-cell-endpoint-fine-d2-x05-state-2026-11-28.json.gz > /tmp/fine-d2-x05-state.json
for m in 1 2; do
  python -W error early_pole_endpoint_continuum_cold_band_pilot.py \
    /tmp/fine-d2-x05-state.json "$m" "/tmp/fine-d2-cold-$m.jsonl"
done
python -W error early_pole_endpoint_cold_band_panel_contrast.py \
  --d1 results/fine-endpoint-continuum-cold-{1,2}.jsonl.gz \
  --d2 results/fine-d2-endpoint-continuum-cold-{1,2}.jsonl.gz \
  --states results/early-pole-cell-endpoint-fine-x05-state-2026-11-28.json.gz \
           results/early-pole-cell-endpoint-fine-d2-x05-state-2026-11-28.json.gz \
  --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-2026-11-28.json
```
