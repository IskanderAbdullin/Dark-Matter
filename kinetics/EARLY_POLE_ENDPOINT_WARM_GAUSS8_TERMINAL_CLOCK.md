# Actual within-G8 terminal pressure time-step check

Recover the G8 1024-step nodal state on x=.04 to .05 from the exact, previously archived G8 x=.04 state. The recovered coarse trace matches the original published trajectory trace exactly on the verified fields. Both clocks then use the same S4 mass for the pressure measurement. Coarse relative to fine changes scalar count by -3.18010e-5 but pressure by only +2.49407e-7 (warm-only pressure: +2.97874e-7).

This is directly measured time-step sensitivity on the last segment from the same initial state. It is much smaller than the matched G8/G7 pressure difference of about 4.99e-4, but is not a global clock error bound and is separate from the conditional count-energy pressure envelope.
The previously recovered mixed255 endpoint gives mixed255/G8 pressure contrasts of +3.3698175e-4 on the two 1024-step states and +3.3697938e-4 on the two 2048-step states, at the same mass. Their difference is only 2.37e-9 in fractional units. Mixed255 and G8 retain their own exact x=.04 initial states in this two-chart comparison.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error early_pole_endpoint_gauss_warm_pressure_clock.py \
  "$p-x05-state-$d.json.gz" "$p-gauss8warm-x04-state-$d.json.gz" \
  "$p-gauss8warm-x05-state-$d.json.gz" "$p-gauss8warm-x05-$d.json.gz" \
  --coarse /tmp/g8-x05-coarse-state.json
```
