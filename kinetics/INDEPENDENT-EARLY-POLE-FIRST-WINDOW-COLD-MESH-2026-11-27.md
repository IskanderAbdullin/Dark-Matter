# Refining cold scalar cells shifts the first-interval gain by less than 0.031%

Use the same joined first interval, hard supports, prepared hard inputs, active unlike response and 16 macro steps as the [time-grid baseline](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-16-REFINEMENT-2026-11-27.md). Change only the number of logarithmic cold scalar cells on the original support [1e−6,1e−5]. Volume overlap sets the initial counts. Coarsening either finer initial mesh recovers the eight-bin counts to below 4e−17 in aggregate absolute difference. The accepted local elastic energy budget and terminal gates are unchanged.

| Cold bins | Native first mean rate | Shift vs 8 | Continuous first mean rate | Shift vs 8 |
| ---: | ---: | ---: | ---: | ---: |
| 8 | .06694243427 | — | .06528878623 | — |
| 16 | .06696254413 | +300.4 ppm | .06527531500 | −206.3 ppm |
| 32 | .06696213698 | +294.3 ppm | .06527666788 | −185.6 ppm |

The 16→32 step changes the gain by about −6.1 ppm and +20.7 ppm in the two inputs. The instantaneous native starting rate is more mesh-sensitive, moving by +.00058577 at 32 vs 8, or .745% of the eight-cell value. These scalar-grid effects do not explain the percent-level difference between the prepared hard inputs on this time grid. This is one scalar support and one time width; it neither establishes the joint time-and-momentum continuum limit nor constrains omitted scalar momenta.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for p in steps cell-conservative; do
  for q in 8 16 32; do
    python -W error early_pole_cold_screened_split_reference.py fine \
      --duration 3.3333333333333335e-7 --count 16 --pair --same-soft-strength 1 \
      --compton-strength 1 --conversion-strength 1 --screening-update full --projection "$p" \
      --snapshots --elastic-energy-budget 5e-15 --scalar-bins "$q" > "$p-q$q.json"
  done
  python -W error early_pole_first_window_cold_bin_refinement_reference.py --projection "$p" \
    "$p-q8.json" "$p-q16.json" "$p-q32.json"
done
```

Raw: `results/early-pole-first-window-cold-mesh-2026-11-27.json.gz`.
