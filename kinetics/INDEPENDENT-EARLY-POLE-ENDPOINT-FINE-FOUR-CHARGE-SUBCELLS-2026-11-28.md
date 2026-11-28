# The fine charge ladder does not yet converge in the late production interval

Keep the 67-bin scalar grid and the [same endpoint-cell collision construction](INDEPENDENT-EARLY-POLE-CELL-ENDPOINT-COLLISION-2026-11-28.md). Split each parent and fermion momentum cell into four subcells before empty-daughter startup. The two- and four-subcell histories both start with the same single 2048-step clock to x=.0075, then retain the finer arm at every restart. The following are the unmodified one-subcell baseline and the independently evolved two- and four-subcell results.

| x | One subcell | Two subcells | Four subcells | Four relative to two |
| --- | ---: | ---: | ---: | ---: |
| .01 | 3.47894426e−7 | 3.47540472e−7 | 3.47431446e−7 | −.0314% |
| .02 | 1.02796002e−5 | 1.02195733e−5 | 1.01987273e−5 | −.2040% |
| .04 | 1.77413501e−4 | 1.73147769e−4 | 1.70375245e−4 | −1.6012% |
| .05 | 3.78441807e−4 | 3.71312973e−4 | 3.60874791e−4 | −2.8112% |

At x=.05 the successive absolute yield decrements have ratio 1.464. They are **larger**, not smaller, under the final doubling, despite the smaller increments seen earlier along the same histories. This rules out treating the observed one-to-two difference as a late spatial error estimate. It does not rule out convergence at finer resolution. The last within-segment time refinement at four subcells changes the final yield by 3.96e−5 relative from a fixed restart; this is a time test of that segment, not a global time error bar. Number, Pauli and independently integrated continuity budgets pass on the reported time samples. Neither agreement with the original archive nor the measured d=1,2,4 sequence warrants an extrapolated limiting yield.

On the common scalar grid, the four-subcell x=.05 histogram has normalized total variation .00666 relative to the original, and the scalar pole mass is .0008646786 (−.846% relative to the original). The new [kinematic threshold interval](results/early-pole-cell-endpoint-fine-d4-threshold-2026-11-28.json) at the upper q edge and cosine zero is 20.526–20.546 in E/T, above the earlier two-subcell interval 18.050–18.114. This extrapolates neither a collision rate nor an actual photon population. The thinner optical interval alone does not remove the large change in the evolved state.

Reproduce from `kinetics/`: run `early_pole_cell_endpoint_startup_pilot.py fine --div 4 --steps 2048 --one` to x=.0075 and save its state; on subsequent calls provide `--restart` and omit `--one`. Use coarse step counts 512, 1024, 1024, 1024 for x=.01, .02, .04, .05 respectively. The exact fine-time result of each segment is archived as `results/early-pole-cell-endpoint-fine-d4-x*-state-2026-11-28.json.gz`. The [continuation receipt](results/early-pole-cell-endpoint-fine-d4-to-x05-2026-11-28.json) preserves both time arms at each restart. Recombine the nested histories using:
```sh
python -W error early_pole_endpoint_charge_ladder_receipt.py \
  results/early-pole-cell-endpoint-fine-{to-x05,d2-to-x05,d4-to-x05}-2026-11-28.json
```
See the [layered comparison](results/early-pole-cell-endpoint-fine-charge-ladder-2026-11-28.json) and [spectrum check](results/early-pole-cell-endpoint-fine-d4-spectrum-2026-11-28.json). The charge ladder and the separate scalar refinement both remain diagnostics of an opt-in benchmark.
