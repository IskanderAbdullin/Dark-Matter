# Match the short hot-case durations before comparing collision channels

The [combined expanding benchmark](FD_PORTAL_COMBINED_ENTROPY_SPLIT_EXPANSION.md) did not compare terminal states for its hot pair because the timestep-selection bounds gave different physical durations. Hold the archived combined initial state and duration at 0.6882487043 GeV and recompute only the timelike counterpart, at both phase-space resolutions. Independently integrate the continuous reference and verify each split trajectory's pressure-work and comoving-entropy budgets.

Combined minus timelike at exactly the same duration is +4.2046e-5 in dark log-temperature, -1.6428e-4 in dark chemical parameter, and +4.3570e-7 in comoving entropy production divided by the initial entropy. The eight-step comparison gives +4.2092e-5, -1.6440e-4 and +4.3578e-7, respectively. Two timelike phase-space meshes agree on the continuous state and associated budgets to within 1.27e-10 in their reported coordinates. The sign of the dark temperature response in this hot test differs from the 0.30 GeV initial state; these are distinct short, deliberately specified thermal states.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
 fd_portal_combined_expansion_matched_horizon.py eos-range.dat \
 results/fd-portal-combined-entropy-split-expansion-2026-11-28.json
```

This numerical test does not predict the full expansion or relic abundance. See the prior tabulated-bath note for EOS input.
