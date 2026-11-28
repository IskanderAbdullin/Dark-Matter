# The G5 versus G4 pressure contrast survives the last-segment clock change

Restart G5 from its own archived (x=.04) nodal state and recompute the .04–.05 segment with 1024 steps. The archived fine 2048-step replay matches all audited trace values exactly. At the common terminal S4-W1 mass, the actual coarse/fine G5 scalar number minus one is −3.18009567e−5 and the corresponding total pressure numerator difference is +2.51663e−7.

| Fixed-mass pressure contrast | Both 1024 steps | Both 2048 steps |
| --- | ---: | ---: |
| G5/G4 minus one | +0.00112598491 | +0.00112598362 |
| G5/G3 minus one | −0.00012930953 | −0.00012932556 |

The change in G5/G4 is 1.29e−9 in absolute relative fraction; that in G5/G3 is 1.60e−8. On this last segment, neither clock change explains the sign or scale of the corresponding warm-chart contrast. This test holds earlier trajectories fixed; it does not bound accumulated clock error, warm-order convergence, or the physical transfer function.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4
w=results/early-pole-geometric-warm-w
d=2026-11-28
python -W error early_pole_endpoint_gauss_warm_pressure_clock.py \
  "$p-x05-state-$d.json.gz" "$p-gauss5warm-x04-state-$d.json.gz" \
  "$p-gauss5warm-x05-state-$d.json.gz" "$p-gauss5warm-x05-$d.json.gz" \
  --coarse /tmp/G5-coarse.json >/tmp/G5-clock.json
gzip -dc "$p-gauss5warm-coarse-x05-state-$d.json.gz" | cmp /tmp/G5-coarse.json -
python -W error early_pole_endpoint_warm_pressure_clock_contrast.py \
  "$w"1-pressure-clock-$d.json "$w"2-pressure-clock-$d.json "$p-gausswarm-pressure-clock-$d.json" \
  --w4 "$w"4-pressure-clock-$d.json --gauss3 "$p-gauss3warm-pressure-clock-$d.json" \
  --gauss4 "$p-gauss4warm-pressure-clock-$d.json" --gauss5 /tmp/G5-clock.json
```
