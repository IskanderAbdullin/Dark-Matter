# Matching the initial global moments does not change the resolved warm gain

Apply the [explicit spectator adjustment](INDEPENDENT-EARLY-POLE-WARM-GLOBAL-MOMENT-CONTROL-2026-11-27.md) as an opt-in initial condition in the short joined q4 distributed-shell control. This is a manufactured comparison. Keep the hard and cold states, scalar coupling and warm quadrature fixed. Use the minimum relative-count change among the 30 remaining eligible spectator markers to restore the original global number, free energy and inverse-energy moment. The adjusted markers stay frozen during the continuation. Their greatest fractional population change is 0.0355%, and the sum of the absolute number changes is (5.46\times10^{-12}).

Before evolution, the original marker gap 0.0008721874799955858 and the evaluated cold rate 0.11773242642703721 are both reproduced at displayed precision. The full initial number mismatch is (2.65\times10^{-23}); terminal energy closure is within (2.23\times10^{-16}) in both arms. The warm gain, (5.229114575239321\times10^{-13}), is unchanged at the reported precision relative to the uncompensated q4 joined test, and remains about 4.87% below the native marker gain. This does not prescribe the physical spectator distribution or imply that a small cold differential is resolved by two separate integrations.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_live_scalar_exterior_decay_reference.py \
  results/early-pole-joined-fine-216-memory-2026-11-27.json.gz \
  --warm-order 4 --warm-compensate --joined --panels 6
```

Raw: `results/early-pole-live-warm-moment-control-2026-11-27.json.gz`.
