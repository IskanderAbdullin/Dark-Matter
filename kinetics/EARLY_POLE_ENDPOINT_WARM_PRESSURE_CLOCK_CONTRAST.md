# Last-segment pressure contrast on all three warm charts

The 1024-step endpoints were recovered from their individual archived .04 restart states, with the same stepper as the selected 2048-step runs. For W1 and W2 all four audited trace columns at eight stored times match their original coarse receipts exactly. The independent warm Gauss clock replay uses the corresponding nodal state and its own two-moment reconstruction. We compare all six measured endpoints at the single W1 x=.05 mass. The cold scalar and charge charts are unchanged.

| Warm representation | Coarse/fine scalar count minus 1 | Coarse/fine pressure numerator minus 1 |
| --- | ---: | ---: |
| W1 geometric | −3.1800977e−5 | +3.03571e−7 |
| W2 geometric | −3.1800959e−5 | +2.83176e−7 |
| Gauss volume | −3.1800955e−5 | +2.54236e−7 |

The difference in the Gauss-vs-W1 pressure contrast between coarse and fine time grids is 4.51663e−8 absolute relative fraction, versus an 8.4489% pressure contrast. For W2-vs-W1 that clock effect is 1.91232e−8 absolute. These gate only the last .04-to-.05 segment holding each chart's previous trajectory fixed, not the cumulative clock error from startup.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_geometric_warm_pressure_clock.py \
  "$p-x05-state-$d.json.gz" "$p-x04-state-$d.json.gz" \
  "$p-x05-state-$d.json.gz" "$p-x05-$d.json.gz" --coarse /tmp/W1-coarse.json
python -W error early_pole_endpoint_geometric_warm_pressure_clock.py \
  "$p-x05-state-$d.json.gz" "$p-warm2-x04-state-$d.json.gz" \
  "$p-warm2-x05-state-$d.json.gz" "$p-warm2-x05-$d.json.gz" \
  --warm --coarse /tmp/W2-coarse.json
python -W error early_pole_endpoint_warm_pressure_clock_contrast.py \
  results/early-pole-geometric-warm-{w1,w2}-pressure-clock-$d.json \
  "$p-gausswarm-pressure-clock-$d.json"
```
