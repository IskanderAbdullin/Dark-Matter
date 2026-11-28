# Warm redshift increment survives halved macro steps and restart

Repeat the [short held-polarization warm-marker experiment](INDEPENDENT-EARLY-POLE-LIVE-WARM-REDSHIFT-2026-11-27.md) at twelve instead of six macro steps, comparing scale held at zero against the same prescribed log(a)=t at exposure 2e−6. The warm gain increases by 3.0016996 and 3.0017129 ppm, respectively, under prescribed redshift. In absolute terms the warm fixed-scale gain changes by about 1.09 ppm on halving the split width; the *matched redshift difference* is much less sensitive. The largest n12 redshift-subtracted total-energy error is below 9.4e−13, still within the unchanged 3e−12 gate. The n6 error is much smaller. These checks do not give a solved cosmological clock or an error bound at larger exposure.

For longer studies the executable now supports atomically restarting inside either matched arm. Checkpoints include the normalized numerical controls, input and archive hashes, solver versions, completed ledger, current state digest and local macro index. A direct two-panel calculation in each arm and a split replay stopping at overall boundaries 1 and 3 (one per arm) produce identical output bytes. No terminal energy gate or step size changes on resume.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
p=early_pole_live_scalar_exterior_frw_reference.py
python -W error "$p" "$a" --panels 12
python -W error "$p" "$a" --panels 12 --expands
python -W error "$p" "$a" --duration 6e-7 --panels 2 --expands > direct.json
python -W error "$p" "$a" --duration 6e-7 --panels 2 --expands --checkpoint native-warm.json --stop-after 1
python -W error "$p" "$a" --duration 6e-7 --panels 2 --expands --checkpoint native-warm.json --resume --stop-after 3
python -W error "$p" "$a" --duration 6e-7 --panels 2 --expands --checkpoint native-warm.json --resume > resumed.json
cmp direct.json resumed.json
```

Raw: `results/early-pole-live-warm-redshift-time-and-restart-2026-11-27.json`.
