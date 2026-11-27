# The evolved H8 cold32 seed already conserves the archived cell numbers

The [frozen cold quadrature comparison](INDEPENDENT-FINITE-MASS-QED-SHARED-POLAR-COLD-MESH-2026-11-27.md)
shows a considerable **hypothetical** scalar weight error if the
native occupation alone is sampled at one node per output cell.
Check whether this construction actually describes the evolved H8
initial condition. It does not.

Recompute every cold32 cell number directly from overlaps of its
q³ interval with the archived native shell intervals, using the
archived signed shell counts with the same nonnegative retention.
Compute the overlap volume independently as a difference of cubes
factored into width times a quadratic. Compare to the initial mode
history of the actual H8 native trajectory.

All retained 32 numbers agree within 8.7e−19 absolute; the total
is .017334279690162456 in both. Replacing this by an occupation
sample at each cell center would instead raise the retained
number by about 7.45% in this fine fixture. That substitution was
not used in the evolved initialization.

For a separate check of the initial inverse-energy weight,
estimate it from the true retained counts at their representative
energies, versus integrating q² f_s(q)/E_s(q) independently on the
native overlap intervals at the same fixed archived scalar mass.
The representative-energy approximation changes the weight by
about −5.23e−8 relative. This last number is only seed initialization
with fixed mass; it is not a bound on errors during coupled
evolution or in an interacting physical rate.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error early_pole_cold32_archived_seed_alias_control.py
```
[Raw actual-seed check](results/early-pole-cold32-archived-seed-alias-control-2026-11-27.json).
