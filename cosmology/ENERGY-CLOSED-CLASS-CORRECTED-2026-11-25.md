# Corrected energy-closed CLASS results — 2026-11-25

After the [gauge pathology control](CLASS-NORMALIZATION-RECOVERED-2026-11-25.md), both synchronous workflows passed the explicit long-wave condition, abundance/early energy checks, independent kinetic grids and momentum sampling/partition tests.

| Production | z=5 half-power k (Mpc^-1) | P/P_CDM at k=1 | P/P_CDM at k=10 |
|---|---:|---:|---:|
| external mother, 480 | .698397 | .318624 | .037106 |
| energy-closed, 480 | .689883 | .311182 | .034978 |

The energy-closed 320→480 half-crossing shifts ~.32%; pre-half mesh drift ~.37%, partition drift <.48%. These are conditional linear model results, not observational constraints. The new [energy-closed run](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45625224293) also tests the hierarchy approximation: at k=.01, default/late/no-fluid ratios agree within .03%; at k=.5 they range .6732–.6775; at k=1 they range .3104–.3600, and at k=1.5 .1907–.2408. Therefore precision of the default fluid-approximation transfer at high k is not established. The quoted k=10 ratio is the default computation only; it must not be assigned sub-percent uncertainty from the momentum-grid tests.

Use the actual hierarchy for convergence around the half crossing and high-k suppression, then repeat for the much earlier and colder x0=.005 case. No new observational likelihood or thermal microphysical consistency is established. The Newtonian-gauge runs from the same day are superseded, not removed from history.
