# Start the cell-endpoint collision solver with empty daughters

The [local snapshot test](INDEPENDENT-EARLY-POLE-CELL-ENDPOINT-COLLISION-2026-11-28.md) does not test the onset of production. Start instead at the original reference x=.005, with the same per-bin parent seed prescription and no scalar or fermion daughters. Fix the node volumes using the charge cell energy moments at this start. The collisional update preserves the same number and Pauli bounds throughout. The evolution retains the original scalar gap and radiation correction, and the energy ledger independently integrates the trace.

At x=.0075 the 1024→2048-step yield changes by 1.94e−4 (fine) and 1.99e−4 (third), down from about .0027 on the first two meshes. The endpoint yields are 5.7796724e−8 and 5.7955002e−8, versus the published archive's 5.7640287e−8 and 5.7856672e−8. Fine's additional 4096-step yield is 5.7799600e−8 (next relative change about 4.98e−5). These are differences between spatial constructions, not all error in either one.

Use the fine 2048-step state at x=.0075 without resetting the populations. The continuation to x=.01 gives a yield of 3.4789443e−7, with 512→1024 time drift 1.44e−7 on that starting state, versus the older archive value 3.4789682e−7. This agreement is conditional on the fixed fine spatial mesh and restart state. It does not establish spatial convergence, the earlier checkpoint's propagated error, or a full evolution to x=.05.

Replay from `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_cell_endpoint_startup_pilot.py fine --steps 1024 --state-path fine-x0075.json
python -W error early_pole_cell_endpoint_startup_pilot.py third --steps 1024 --state-path third-x0075.json
python -W error early_pole_cell_endpoint_startup_pilot.py fine --steps 4096 --one
python -W error early_pole_cell_endpoint_startup_pilot.py fine --x-end .01 --steps 512 --restart fine-x0075.json
```
The [combined receipt](results/early-pole-cell-endpoint-startup-replay-2026-11-28.json) records both time ladders and the restart. The compressed state files contain node occupations and their grid identifiers; decompress to JSON before passing a state to the replay. No production driver is replaced.
