# Direct continuation through the local slow-mode clock region

Continue the exact tabulated-EOS direct timelike and combined states from log(a/a_start)=0.25 to 0.55. Keep all three physical ledgers normalized to the original 0.90 GeV state, rather than renormalizing at the restart. Independent BDF/Radau and phase-space meshes match the reported coordinates within 1.84e-8; the maximum normalized conservation defect is 5.99e-10. The terminal visible temperature is about 0.53643 GeV. No steady-state assumption is used to evolve through the temperatures where previous local rate ratios crossed unity.

At the terminal clock, combined minus timelike is +3.5932e-4 in dark chemical parameter, -4.2540e-5 in dark log-temperature, and +2.16045e-4 in comoving number divided by the original number. The combined-minus-timelike accumulated entropy difference changes sign between the sampled 0.25 and 0.30 expansion clocks. At 0.55 it is +2.98693e-8 of original entropy. Both cumulative entropy productions remain positive, so the sign change concerns a comparison across distinct solutions, not the sign of a physical source.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_direct_expansion_tracking_continuation.py \
 eos-range.dat results/fd-portal-direct-combined-expansion-tracking-2026-11-28.json
```

Use the same local EOS input as the prior tabulated-bath benchmark. This is a finite direct numerical continuation with restricted vacuum charged-lepton rates, not a complete cosmological history or a decoupling and relic abundance determination.
