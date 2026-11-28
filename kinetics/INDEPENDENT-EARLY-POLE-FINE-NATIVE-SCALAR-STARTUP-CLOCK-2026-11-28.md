# The cold shape discrepancy survives matching the startup time resolution

The initial [scalar resolution comparison](INDEPENDENT-EARLY-POLE-FINE-COLD-SCALAR-RESOLUTION-2026-11-28.md) used the existing 2048-step unrefined startup alongside new refined histories with a 4096-step startup. To isolate that difference, rerun the unrefined 67-shell scalar chart at the same fine two-subcell charged resolution using a single 4096-step empty-daughter startup to x=.0075. Then retain each later segment's original time mesh counts and propagate the new state through all four joins to .05. Both unrefined histories have the same scalar and charged charts and reaction parameters.

| x | Relative yield change from startup refinement | TV of normalized unrefined scalar spectrum |
| --- | ---: | ---: |
| .0075 | +4.975e−5 | 6.228e−5 |
| .01 | +9.635e−6 | 6.235e−5 |
| .02 | +7.264e−7 | 6.217e−5 |
| .04 | +4.592e−7 | 6.204e−5 |
| .05 | +9.183e−7 | 6.191e−5 |

At x=.05, all three scalar charts can now be compared using identical 4096-step starts and identical later step counts. On common native bins, the unrefined→S2 comparison still gives **8.015%** normalized TV, versus **0.00619%** from doubling the unrefined startup clock alone. The difference in startup resolution therefore does not account for the observed coarse spectral difference. This comparison does not bound the full trajectory's time error or turn any one spectrum into a converged physical prediction. In particular it does not measure the independent phase error of the last two refined spectra.

Reproduce using `early_pole_cell_endpoint_startup_pilot.py fine --div 2 --steps 4096 --one` and saving the x=.0075 endpoint. From each exact restart use coarse counts 512, 1024, 1024, 1024 to x=.01, .02, .04, .05. This replay used the [verified optional exact capacity path](INDEPENDENT-EARLY-POLE-ENDPOINT-CAPACITY-FAST-PATH-2026-11-28.md); the reference path has the same update rule. The [continuation receipt](results/early-pole-fine-charge2-native-scalar-start4096-continuation-2026-11-28.json), [clock shape audit](results/early-pole-fine-charge2-native-scalar-start-time-shape-2026-11-28.json), and [fully matched x=.05 comparison](results/early-pole-fine-charge2-scalar-matched-clock-x05-2026-11-28.json) archive the results and both exact state sequences are available alongside them.
