# Locate the crossing of two accumulated entropy histories

Use the exact two archived physical states at log(a/a_start)=0.25 and independently integrate both through 0.30 with dense output. Inspect 41 sample clocks to isolate the measured sign reversal, bracket the root, then verify the difference's positive derivative there. The two phase-space orders and the BDF/Radau controls locate the sampled crossing at log(a/a_start)=0.2902854, with an absolute clock difference of 1.12e-8. Both visible temperatures are about 0.68403 GeV at that point.

Only the **difference** between two accumulated entropy histories crosses zero. Neither cumulative production becomes negative, and there is no identified physical threshold. The precision applies to the two numerical models and the declared initial states, not to an early-universe prediction.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_accumulated_entropy_history_crossing.py \
 eos-range.dat results/fd-portal-direct-combined-expansion-tracking-2026-11-28.json \
 results/fd-portal-direct-expansion-tracking-continuation-2026-11-28.json
```

Use the prior local authors' EOS normalization; no copy of the EOS table is included.
