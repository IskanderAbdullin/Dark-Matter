# Third-grid endpoint-cell evolution through x=.05

The [fine continuation](INDEPENDENT-EARLY-POLE-CELL-ENDPOINT-FINE-CONTINUATION-2026-11-28.md) is now paired with the third-grid replay. Start from the saved empty-daughter solution at x=.0075, and use the entire node-occupation state at each join through x=.01, .02, .04 and .05. Each segment has two time meshes at the same fixed initial state and unchanged spatial representation.

| x | Third endpoint-cell scalar yield | Relative to original third archive | Last fixed-start time change |
| --- | ---: | ---: | ---: |
| .01 | 3.4741265e−7 | −.001192 | 1.36e−7 |
| .02 | 1.0257438e−5 | −.005226 | 1.74e−6 |
| .04 | 1.7490652e−4 | +.021595 | 9.55e−5 |
| .05 | 3.7652201e−4 | +.034138 | 2.47e−5 |

The final fine and third endpoint-cell yields therefore differ by about .51%, down from 1.67% between the two original results. This is a comparison of two complete fixed-grid trajectories, not a spatial convergence bound. The original and new distributions on the third scalar grid have normalized total variation .00601; the new scalar pole rises by 1.125%. The cold abundance remains strongly concentrated below q=1e−5, but final spectrum shape and downstream reaction rates are not certified. Pauli bounds and the independently integrated energy and number ledgers passed on the reported time samples.

Replay a segment by decompressing the prior `early-pole-cell-endpoint-third-*-state-2026-11-28.json.gz` into JSON. For example, from `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_cell_endpoint_startup_pilot.py third \
 --x-end .05 --steps 1024 --restart WORK_X04_STATE.json --state-path WORK_X05_STATE.json
python -W error early_pole_endpoint_final_spectrum_audit.py WORK_X05_STATE.json
```
The [continuation receipt](results/early-pole-cell-endpoint-third-to-x05-2026-11-28.json) and [spectrum comparison](results/early-pole-cell-endpoint-third-spectrum-x05-2026-11-28.json) record the fixed-chart results. No default driver or cosmological prediction is changed.
