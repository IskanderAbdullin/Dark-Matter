# The G6/G5 pressure contrast survives the final-step clock change

The coarse G6 x=.05 state is independently recomputed from the hashed x=.04 restart of the full G6 chain. This holds the previous G6 trajectory fixed and uses 1024 steps instead of the archived 2048; all recorded coarse trace columns replay exactly. The coarse/fine changes in G6 number and pressure at the shared native S4 measurement mass are −3.18010e−5 and +2.52860e−7, in relative fraction.

The measured G6/G5 pressure ratio minus one is −0.000817581633 on the paired 1024-step endpoints and −0.000817582829 on the paired 2048-step endpoints. The absolute fractional difference is 1.1961e−9. Together with the separate first-continuation control, this tests two distinct intervals, not the remaining clock stages, the total time discretization error, or a physical pressure prediction.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-geometric-warm-w
python -W error early_pole_endpoint_gauss_warm_pressure_clock.py \
 "$p-x05-state-$d.json.gz" "$p-gauss6warm-x04-state-$d.json.gz" \
 "$p-gauss6warm-x05-state-$d.json.gz" "$p-gauss6warm-x05-$d.json.gz" \
 --coarse /tmp/G6-coarse.json >/tmp/G6-clock.json
python -W error early_pole_endpoint_warm_pressure_clock_contrast.py \
 "$w"1-pressure-clock-$d.json "$w"2-pressure-clock-$d.json "$p-gausswarm-pressure-clock-$d.json" \
 --w4 "$w"4-pressure-clock-$d.json --gauss3 "$p-gauss3warm-pressure-clock-$d.json" \
 --gauss4 "$p-gauss4warm-pressure-clock-$d.json" --gauss5 "$p-gauss5warm-pressure-clock-$d.json" \
 --gauss6 /tmp/G6-clock.json
```
