# Temperature equality with nonzero chemical force in the direct FD path

The independent dense direct FD histories both reach equal dark and visible temperatures around 0.434 GeV. In the combined path this occurs at log(a/a_start)=0.780384, with chemical parameter eta ≈ -0.06940. The elastic heat vanishes at exact common temperature, but the timelike number and heat terms remain nonzero and positive total entropy is produced. The branches use two phase-space orders and BDF/Radau controls; the larger independent difference in the crossing coordinates is 2.60e-8. Equal bath temperatures are not equivalent to chemical equilibrium.

Separately, compare the equilibrium linearized **instantaneous source** to the exact FD source on the four already verified local-clock temperature states. At the combined slow label the linearized number, heat and entropy magnitudes are higher by 0.775%, 0.608% and 1.30%; at the combined fast label they are higher by 13.20%, 15.53% and 14.22% respectively. All four cases have the same source signs before and after linearization. These source errors are distinct from the much larger error of imposing a stationary lag along the trajectory.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28
python -W error fd_portal_direct_temperature_equality_with_chemical_lag.py \
 eos-range.dat "results/fd-portal-direct-expansion-055-096-$d.json"
python -W error fd_portal_live_clock_source_linearity.py \
 "results/fd-portal-live-states-at-local-clock-temperatures-$d.json"
```

Use the earlier local EOS normalization. These are properties of finite, restricted direct FD trajectories and selected local comparisons, not thermal or chemical decoupling determinations.
