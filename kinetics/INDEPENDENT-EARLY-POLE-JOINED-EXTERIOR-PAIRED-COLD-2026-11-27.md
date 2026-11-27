# On a shared grid the joined warm tail has a very small cold feedback

Using the same fine native final n216 snapshot and warm marker selection q∈[.1,5] as the [joined warm continuation](INDEPENDENT-EARLY-POLE-EXTERIOR-JOINED-AND-DISTRIBUTED-2026-11-27.md), evolve a baseline and the explicit difference from allowing warm markers to evolve. The total state differences are not estimated by subtracting separately integrated cold populations. Slow rates are stepped on one RK4 grid, differencing each common hard channel before summing its response. Each backward elastic update refreshes the active unlike response in its own arm. The hard difference is updated from the two elastic increments, while the elastic solve does not change the scalar counts.

| Joined macro steps | Baseline cold gain | Added warm gain | Warm feedback on cold |
| --- | ---: | ---: | ---: |
| 6 | 2.3563420188e−7 | 5.4968656008e−13 | −7.48968e−19 |
| 12 | 2.3563759048e−7 | 5.4968715599e−13 | −7.49155e−19 |

The explicit differential changes by about .025% on halving the time step and is about −3.18e−12 of the baseline cold gain. The warm gain matches the separately advanced joined-marker integrator within 1.2e−23 absolute at each resolution. The terminal hard-count and total-energy defects stay within 3e−12 in both arms. At n12 their energy drifts are about +9.22e−13 and are nearly common to both arms. These ledgers and the split comparison do not establish a general error bound for the tiny differential, nor predict an abundance, a full warm reconstruction, or a long-time result.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
for n in 6 12; do
  python -W error early_pole_live_scalar_exterior_paired_difference_reference.py "$a" \
    --joined --duration 2e-6 --steps "$n" --order 8
done
```

Raw: `results/early-pole-joined-exterior-paired-cold-2026-11-27.json`.
