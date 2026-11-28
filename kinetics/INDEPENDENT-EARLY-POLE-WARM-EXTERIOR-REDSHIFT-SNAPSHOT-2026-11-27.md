# Reciprocal warm packets on matched prescribed redshift snapshots

Prepare a prerequisite for a live expanding warm tail. Keep the final fine n216 hard and scalar comoving populations fixed. Evaluate them at imposed comoving scale log(a) = 0, .00002, .000072, .0003 and .001, using the checked hard redshift and spectator-trace machinery from the [joined FRW study](INDEPENDENT-EARLY-POLE-FRW-LONG-EXPANDING-CLOCKS-2026-11-27.md). The reciprocal band kernel now optionally uses the same a-scaled bare mass in its vacuum decay prefactor as the established cold-window kernel, with the default a=1 unchanged. The charged occupation maps and scalar gap calculation otherwise stay as in their respective fixed-snapshot state.

Across the two hard reconstructions and all five scales, the maximum absolute difference from the independent original cold-window mode rates is 1.05e−12 or less. The warm exchange has absolute reciprocal energy and charge-number defects below 2e−22 at every scale and positive entropy production. At log(a)=7.2e−5, the q=.1–5 marker warm rate is about 216.14 ppm higher than at a=1 on the held comoving state. Almost all of this snapshot change comes from the explicit vacuum prefactor (216.06 ppm at the same shifted hard energy), rather than the change of the hard energy map. This does not integrate a cosmological clock, evolve the populations, or account for a redshift-dependent medium response.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for p in steps cell-conservative; do
  python -W error early_pole_exterior_frw_snapshot_reference.py \
    results/early-pole-joined-fine-216-memory-2026-11-27.json.gz --projection "$p"
done
```

Raw: `results/early-pole-warm-exterior-redshift-snapshot-2026-11-27.json`.
