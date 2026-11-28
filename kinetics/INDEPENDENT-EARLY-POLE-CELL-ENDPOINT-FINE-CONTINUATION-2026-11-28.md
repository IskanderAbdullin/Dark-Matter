# The fine endpoint-cell evolution reaches the original final scale

Continue the [fine early replay](INDEPENDENT-EARLY-POLE-CELL-ENDPOINT-STARTUP-2026-11-28.md) through x=.02, .04 and .05. At every join use the complete preceding saved state, with no changes to endpoint-node volumes, grid dimensions or occupations. Independently halve the time step on each segment while holding its starting state fixed. The production solver in the original archive is not replaced.

| x | Fine endpoint-cell scalar yield | Relative to original fine archive | Last local time change |
| --- | ---: | ---: | ---: |
| .01 | 3.4789443e−7 | −6.87e−6 | 1.43e−7 |
| .02 | 1.0279600e−5 | −.003910 | 1.23e−6 |
| .04 | 1.7741350e−4 | +.010734 | 6.48e−5 |
| .05 | 3.7844181e−4 | +.022289 | 1.61e−5 |

The fixed-start time refinements do not bound error propagated through earlier restart states. The independently summed energy and charge budgets pass on each segment and Pauli bounds hold; the final maximum charged occupation is below .499. The change from the previously archived run is a difference between spatial/kinetic representations, not merely the measured temporal drift.

At x=.05 both snapshots use the same scalar edges. Directly summing shell counts gives 1.78% total variation between their *normalized* scalar histograms and a +.7375% shift of the scalar pole mass. Both retain over 99.998% of the number below q=1e−5, but that coarse abundance statement does not constrain a small warm tail or the cold reaction aperture. The new snapshot is not spatially converged and is not a CLASS input.

To replay a segment, first decompress the prior `early-pole-cell-endpoint-fine-*-state-2026-11-28.json.gz` into a work path. From `kinetics/`, run:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_cell_endpoint_startup_pilot.py fine \
  --x-end .05 --steps 1024 --restart WORK_X04_STATE.json --state-path WORK_X05_STATE.json
python -W error early_pole_endpoint_final_spectrum_audit.py WORK_X05_STATE.json
```
The [continuation receipt](results/early-pole-cell-endpoint-fine-to-x05-2026-11-28.json) preserves all matched segment reports. The [spectrum comparison](results/early-pole-cell-endpoint-fine-spectrum-x05-2026-11-28.json) is calculated from the complete archived final occupations.
