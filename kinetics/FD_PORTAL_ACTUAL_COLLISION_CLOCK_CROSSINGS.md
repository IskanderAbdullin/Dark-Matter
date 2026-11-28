# Unity crossings of the live fixed-energy FD collision rates

The earlier reference clock temperatures came from near-equilibrium collision rates at specified visible temperatures. Here locate where the instantaneous collision Jacobian, evaluated on the exact directly evolved finite state at fixed total energy, has an eigenvalue of magnitude H. The root is along the expanding history, but the Jacobian differentiated here is collision-only. At each trial point recompute the finite-state slopes in dark number and energy; do not simply interpolate between the previously reported rate values.

| Restricted channel | Mode | x = log(a/a_start) | Visible T (GeV) | T above earlier reference (MeV) |
|---|---|---:|---:|---:|
| Timelike | Slow | 0.50954303 | 0.55699929 | 0.8511 |
| Timelike | Fast | 0.89643955 | 0.39098857 | 5.8064 |
| Timelike plus elastic scattering | Slow | 0.51754070 | 0.55286570 | 0.9897 |
| Timelike plus elastic scattering | Fast | 0.90352984 | 0.38849660 | 6.0873 |

All four occur earlier and at higher visible temperatures than their original equilibrium-based reference labels. In the combined case the fast label shifts by -0.01754489 in log(a/a_start), about an order of magnitude larger than the slow shift (-0.00192480). Two independent BDF/Radau dense trajectories and phase-space/finite-difference resolutions agree within 5.92e-10 in the most discrepant reported coordinate. The rate signs at five samples bracket each crossing in the listed intervals (.50–.55 or .85–.91). Entropy, expansion work and dark number ledgers retain their previous gates.

Comparing the archived fast-root state with the independently extended final state at x=1.20, the comoving count rises by another 1.13786% (timelike) or 1.10299% (combined) of the original count. Evaluate the same archived FD density as exp(3x) n(T_dark, eta) / n(0.90 GeV, 0) at each endpoint; this does not require interpolating the source after the root.

These are isolated, directly bracketed local rate/H roots for the specified vacuum charged-lepton test, not proof of a global unique mode crossing, a freeze-out temperature, or a relic-abundance prediction.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fd-portal; e=eos-range.dat
python -W error fd_portal_actual_collision_clock_crossing.py "$e" \
 "$p-direct-expansion-tracking-continuation-$d.json" \
 "$p-direct-expansion-055-096-$d.json" \
 "$p-live-states-at-local-clock-temperatures-$d.json"
```

Use the exact normalized public EOS in `FD_PORTAL_TABULATED_BATH_ENTROPY_AND_MODES.md`.
