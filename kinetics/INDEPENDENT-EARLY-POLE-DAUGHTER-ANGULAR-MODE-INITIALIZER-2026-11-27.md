# A positive discrete scalar initializer on the shared angular cuts

The [shared angular test](INDEPENDENT-EARLY-POLE-DAUGHTER-COLD-VOLUME-LADDER-2026-11-27.md)
evaluates a static source. For a possible separately tested
trajectory, represent the 39 common scalar volume intervals
with two independent positive-weight nodes each. Intersect the
panels with the original scalar density jumps and initialize
each node's number as its positive quadrature weight times
the native phase-volume density. The archived spectator count
array does not change. Reinitialize the same scalar pole with
the existing discrete gap prescription; this is not an independent
measurement of a scalar mass or coupling.

The two hard fixtures give exactly equal node positions,
weights and seeds in all 78 modes. Regrouped counts match
the original eight-bin seed within 1.74e-18. The evaluated
pole displacement from the old seed is 2.06e-18. With the
*reconstructed* hard occupations and each arm's own dressed
geometry, the local-minus-native initialized rate is
-0.01629912710 at hard depth 1 and -0.01227352390 at depth 4.
These are not the rates for the distinct fixed-geometry
pointwise archive comparison.

Reproduce with `early_pole_daughter_angular_mode_initializer.py`
on the original and depth-4 eight-mode fixtures and
`results/early-pole-daughter-shared-cold-angular-sampling-2026-11-27.json`.
This is an initialized positive representation only. It does
not show that the modes remain admissible or that a later
cold gain agrees on this mesh under an actual time solver.
