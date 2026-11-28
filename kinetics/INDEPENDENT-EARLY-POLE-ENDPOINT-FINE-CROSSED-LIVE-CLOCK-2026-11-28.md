# Similar diagonal decay clocks conceal large opposing state changes

The [full fine charge ladder](INDEPENDENT-EARLY-POLE-ENDPOINT-FINE-FOUR-CHARGE-SUBCELLS-2026-11-28.md) changes both the charged state and, indirectly, the scalar occupation. All three continue on the same 67 scalar bins. At each of x=.04 and x=.05, keep the evolved charged occupations on their own chart and replace only the scalar occupation array by that of each of the three evolved cards in turn. No charged remap is required and the scalar geometry is exactly common. Recompute the scalar mass and expansion-clock factor for every recombined state. Compare the instantaneous logarithmic number source, `H*sum(scalar changes)/scalar number`, including both stimulated gain and inverse loss. This is from the evolution's decay kernel, not the previously measured fixed-momentum transverse photon coefficient.

At x=.05 the complete crossed table is:

| Charge chart and evolved charges ↓ / evolved scalar state → | d=1 | d=2 | d=4 |
| --- | ---: | ---: | ---: |
| d=1 | 3.68670 | 7.26257 | 8.26814 |
| d=2 | .48741 | 3.36577 | 7.12404 |
| d=4 | .14676 | .19356 | 3.24607 |

The d=2 to d=4 diagonal change is −.11970. Holding the d=2 scalar state while changing the charge chart **and** its evolved occupations gives −3.17221. Holding the d=2 charge chart/state while changing only the scalar snapshot gives +3.75827. The finite interaction (subtract both single replacements from the diagonal change) gives −.70576. At x=.04 the corresponding terms are −2.21156, +2.24720 and −.73238, summing to −.69674. This difference-of-differences depends on the chosen base state and does not identify a unique causal contribution to the full history.

Do not interpret small changes along the diagonal as evidence that either sector is already stable. The mixed cards are frozen counterfactuals, not physical states attained by replay. Event energy and pair counts close; the relative quadrature drift between radial orders 4 and 8 is at most 4.22e−9 on the reported gross and net results. The scope excludes rates integrated over cosmological time or portal processes outside the diagnostic decay.

See the [x=.04](results/early-pole-endpoint-fine-crossed-live-clock-x04-2026-11-28.json) and [x=.05](results/early-pole-endpoint-fine-crossed-live-clock-x05-2026-11-28.json) receipts, which include gross forward/inverse terms and the fitted mass and clock factor. From `kinetics/`, run with the decompressed fine d=1,2,4 states at the same x:
```sh
OPENBLAS_NUM_THREADS=1 python -W error early_pole_endpoint_crossed_live_clock.py \
  FINE_D1_STATE.json FINE_D2_STATE.json FINE_D4_STATE.json
```
