# Propagate improved time steps through every fine restart

The [first complete fine chain](INDEPENDENT-EARLY-POLE-CELL-ENDPOINT-FINE-CONTINUATION-2026-11-28.md) refined each segment at a fixed incoming checkpoint. Here halve the previously retained fine step on **every** segment, beginning with the independently run 4096-step startup. Use only the newly obtained state at the next restart; never switch back to a coarser state. The charge and scalar discretizations remain unchanged.

| Checkpoint x | Refined scalar yield | Relative change from prior fine chain |
| --- | ---: | ---: |
| .01 | 3.4789786e−7 | +9.87e−6 |
| .02 | 1.0279611e−5 | +1.08e−6 |
| .04 | 1.7741641e−4 | +1.64e−5 |
| .05 | 3.7844527e−4 | +9.15e−6 |

Both chains hold the Pauli, number and independently integrated energy conditions. This specific global time refinement is much smaller than the spatial change [on doubling charge intervals](INDEPENDENT-EARLY-POLE-ENDPOINT-FINE-CHARGE-REFINEMENT-2026-11-28.md), −1.884% at x=.05, so the measured spatial sensitivity cannot be attributed to the tested restart/time error. This does not establish a general time-error bound or spatial convergence.

Regenerate the refined chain using `early_pole_cell_endpoint_startup_pilot.py fine --one --steps 4096` for the startup, `--x-end .01 --one --steps 2048` on its saved state, and then `--one --steps 4096` on each subsequent segment to .02, .04 and .05, always continuing the preceding state. From `kinetics/`, pass the four raw JSON reports to `early_pole_endpoint_chain_time_audit.py` after `results/early-pole-cell-endpoint-fine-to-x05-2026-11-28.json`. The [receipt](results/early-pole-endpoint-whole-chain-time-audit-2026-11-28.json) includes all fine-time reports; the final exact occupation state is also archived.
