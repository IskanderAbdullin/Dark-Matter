# The simultaneous held-screening result is stable under a tighter solver tolerance

Repeat the [simultaneous charged-logit calculation](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-UNSPLIT-CHART-2026-11-27.md) at the same endpoints, hard and cold meshes, fixed screening and model inputs. Reduce the DOP853 relative tolerance from 3e−9 to 3e−10 and the baseline absolute tolerance from 2e−14 to 2e−15. Both calculations keep the same report grid and limit on the maximum step width.

| Prepared input | Nominal first mean | Tighter first mean | Relative change |
| --- | ---: | ---: | ---: |
| Native steps | .0664989389804 | .0664989389595 | −3.13e−10 |
| Continuous hard reconstruction | .0655234826188 | .0655234826501 | +4.77e−10 |

The tighter relative contrast is 1.4668750% of the native mean. Across the common report samples the maximum absolute count differences are 1.04e−17 and 1.39e−17. Terminal energy, entropy and hard count pass the unchanged gates. The extra integration calls test solver tolerance only; this does not reduce phase-space discretization, frozen-screening or preparation error, and does not validate a physical cosmology.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
for p in steps cell-conservative; do
  python -W error early_pole_joined_simultaneous_chart_reference.py "$a" \
    --projection "$p" --rtol 3e-10 --base-atol 2e-15
done
```

Raw: `results/early-pole-first-window-unsplit-tolerances-2026-11-27.json`.
