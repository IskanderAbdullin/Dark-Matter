# Skip the capacity limiter only when it cannot alter the step

The [bounded endpoint-cell update](INDEPENDENT-EARLY-POLE-CELL-ENDPOINT-COLLISION-2026-11-28.md) limits signed reaction extents before applying them, rather than clipping individual occupations after deposition. The same arithmetic can take a shorter path in the frequent case when **no direction is limited**. Compute the identical gross forward/inverse arrays, signed microscopic extents, and outgoing/incoming capacity ledgers. If every outgoing demand is at most the old safety margin times that cell's inventory, and every incoming demand is at most the same margin times that cell's hole capacity, then every finite ratio used by the original limiter is ≥1. Zero demand gets the original default multiplier of one. All original per-event factors therefore clip to one, including the unbounded scalar incoming direction.

In this case, accumulate the same signed arrays without multiplying every term by one, and omit the unused first net-source scatter. If either sufficient condition fails, delegate the whole step to the original limiter with the same rate state and time increment. Never skip limiting on the basis of typical amplitudes or a previous step. The default startup path still uses the original implementation; the optional `--fast-capacity` switch changes neither the time method, the charge grid, nor the physics.

The [fixture](results/early-pole-endpoint-capacity-fast-fixture-2026-11-28.json) compares small steps, strongly active limiter stress cases and short evolving histories at matched resolution. The resulting occupations and reported limiter/budget data are bitwise identical on these tests. Its CPU timing compares the fast and original implementations on matched inputs in the same process. These timings are illustrative and do not measure the acceleration of a full production trajectory. A step that actually needs limiting may be slower because the original limiter still runs after the sufficient-condition test.

From `kinetics/`:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error early_pole_endpoint_capacity_fast_fixture.py
python early_pole_cell_endpoint_startup_pilot.py --help
```
