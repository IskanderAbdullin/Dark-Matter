# Which native mesh axis moves the archived decay clock?

The [cross-mesh endpoint comparison](INDEPENDENT-EARLY-POLE-ENDPOINT-STENCIL-REMAP-2026-11-28.md) combines changes across three different mesh axes. On each frozen input independently, vary the destination scalar, fermion daughter and charged parent grid sizes separately. Set each axis to its fine or third resolution, conservatively remap each species in p³, and evaluate the same instantaneous original 1↔2 driver at x=.05. There are eight target meshes per archive.

For the log growth clock, average the marginal contribution from changing each axis from third to fine over all six possible orderings. The log contributions sum to the full native target-axis change, and their exponentials are the multiplicative factors below.

| Frozen input | Scalar axis | Fermion daughter axis | Parent axis |
| --- | ---: | ---: | ---: |
| Fine | 1.00540 | 1.50976 | 1.02942 |
| Third | 1.01270 | 2.22600 | 1.04840 |

The daughter representation accounts for the largest sensitivity in both fixed-state experiments. This points to a targeted daughter-cell refinement at fixed scalar and parent meshes before spending effort on much tighter time integration. The result combines the target discretization and the finite-cell re-representation of the same input. It does not show that either historical solution is converged or that the daughter grid alone explains the difference in their evolving total counts.

Replay `OPENBLAS_NUM_THREADS=1 python -W error early_pole_endpoint_axis_design_pilot.py` from `kinetics/` and compare the [receipt](results/early-pole-endpoint-axis-design-2026-11-28.json).
