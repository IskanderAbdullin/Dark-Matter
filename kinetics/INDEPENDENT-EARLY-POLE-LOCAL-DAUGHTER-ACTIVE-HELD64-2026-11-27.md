# Active response refresh keeps the positive 64-step local effect

At the original fine supports, reproduce the local prepared case
through 64 split steps with active response refresh. Compare to
the matched native path with the same response rule, and to
the two equally resolved held-response runs. Keep initial hard
and cold data, exposure, collision switches, wave geometry and
energy budgets fixed. Both new paired paths pass global closure
and the local elastic gate (maximum individual defect 4.45e-15
under 5e-15).

| Outer time intervals | Held local minus native | Active local minus native |
|---|---:|---:|
| 16 | -1.00295605e-5 | -1.00299144e-5 |
| 64 | +3.84570917e-4 | +3.84571749e-4 |

At 64 steps the positive local effect persists with active
screening. The active-minus-held change is +8.33e-10 in mean
rate (about 2.17 ppm of the held shift); the largest absolute
change in paired cumulative count over the report grid is
8.09e-16. Both 64-step curves cross in the same saved interval,
1.667e-7 to 1.875e-7. Thus the sign change with split refinement is
not explained by whether the implemented response is held or
refreshed at these widths. This does not establish the time limit
of the active model, equivalence of absolute active and held
rates, or the effect of a different physical dielectric response.

For the local active control use
`early_pole_prepared_hard_swap_trajectory_reference.py` with
`--case local-daughter --steps 64` (no `--hold-screening`). The
four matched widths and simultaneous reference are checked by
`early_pole_local_daughter_time_control_reference.py`, with its
optional `--active64` native/local pair.

The active runs actually refresh on all 64 steps, versus zero
refresh events in the held runs. The native and local final mean
soft-packet ratios are 0.999718598 and 0.999719345. The absolute
mean cold rates change by -3.894e-8 and -3.810e-8, respectively,
so most of the change cancels in the local-minus-native contrast.
The maximum recorded response changes are 0.0849 and 0.0743 of
the initial squared Debye mass; that pointwise diagnostic is not
the same as the average packet-weight change.
