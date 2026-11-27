# Cold momentum and time refinements do not add exactly

Use the same prepared first-window model and original cold interval. Cross 16 versus 64 time steps with 8 versus 32 cold bins. The 16-step spatial and eight-bin temporal arms come from the [cold mesh](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-COLD-MESH-2026-11-27.md) and [time-local memory](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-TIME-LOCAL-MEMORY-2026-11-27.md) controls. The 64-step, 32-bin trajectories are new. Separate tests verify identical initial states across widths within each momentum grid, nested initial counts across momentum grids, and unchanged terminal conservation and step-budget gates.

| Time steps / cold bins | Native first mean | Continuous first mean | Contrast / native |
| --- | ---: | ---: | ---: |
| 16 / 8 | .06694243427 | .06528878623 | 2.4703% |
| 64 / 8 | .06661149218 | .06546547131 | 1.7205% |
| 16 / 32 | .06696213698 | .06527666788 | 2.5170% |
| 64 / 32 | .06661735825 | .06545627236 | 1.7429% |

At 16 time steps, changing the cold mesh gives +294 and −186 ppm of the corresponding gain; at 64 it gives +88 and −141 ppm. The signed time–momentum cross difference, relative to the coarsest gain, is −207 ppm (native) and +45 ppm (continuous). Thus one should not add the coarse time and cold-bin corrections as if they were independent. Within these tested grids the prepared-input effect remains substantially larger than the cold mesh correction. This finite cross experiment does not estimate a joint continuum limit or address the omitted infrared or warm support.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for p in steps cell-conservative; do
  for n in 16 64; do
    for q in 8 32; do
      python -W error early_pole_cold_screened_split_reference.py fine \
        --duration 3.3333333333333335e-7 --count "$n" --pair --same-soft-strength 1 --compton-strength 1 \
        --conversion-strength 1 --screening-update full --projection "$p" --snapshots \
        --elastic-energy-budget 5e-15 --scalar-bins "$q" > "$p-n$n-q$q.json"
    done
  done
  python -W error early_pole_joined_time_and_cold_mesh_reference.py --projection "$p" \
    "$p-n16-q8.json" "$p-n64-q8.json" "$p-n16-q32.json" "$p-n64-q32.json"
done
```

Raw: `results/early-pole-first-window-time-cold-mesh-cross-2026-11-27.json.gz`.
