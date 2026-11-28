# A resolved shape gap need not identify the direction of a mean momentum change

Hold each chart's exact gross mass in each native cold panel fixed, and allow arbitrary subpanel location. For a monotone observable, summing its left-endpoint and right-endpoint values weighted by the native masses gives sharp lower and upper expectations, approached by concentrating inside the endpoints. This also bounds the RMS momentum using the monotone square and square root. For speed use a single positive numerical mass, the fine D2 scalar gap, on all grids. These are response-weighted moments, not particle number moments.

Across all four complete chart pairs, the panel ledgers permit both signs in right-minus-left mean momentum, RMS momentum, and mean speed for each of forward number, inverse number, and forward excess. In the forward channel:

| Pair | Right minus left mean q | Right minus left mean speed at common mass |
| --- | ---: | ---: |
| fine D2 / third D2 | [−6.483e−7,+3.307e−7] | [−7.426e−4,+3.788e−4] |
| S2 / S4 | [−3.551e−7,+3.074e−7] | [−4.068e−4,+3.521e−4] |
| S4 / S8 | [−2.064e−7,+1.963e−7] | [−2.364e−4,+2.249e−4] |
| S8 / S16 | [−1.156e−7,+1.136e−7] | [−1.324e−4,+1.302e−4] |

Both quadratures agree to 4.08e−17 or better in the difference interval endpoints. Even the resolved nonnested fine/third shape gap does not determine the sign of the mean shift from fixed panel masses alone. Additional within-panel information would be needed before drawing an inference about relative momentum or velocity. No physical transfer or evolution is inferred from this diagnostic.

```sh
cd kinetics
python -W error early_pole_endpoint_cold_panel_moment_bounds.py --root results \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-2026-11-28.json
```
