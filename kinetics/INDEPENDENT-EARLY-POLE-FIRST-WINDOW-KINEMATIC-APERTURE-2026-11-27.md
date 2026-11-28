# The cold packet geometry itself selects the same two hard cells

Independently find the parent-momentum intervals that meet the energy and angular constraints at the saved 64-step first-window endpoints. This calculation uses the variational cell energy maps and actual scalar pole but no occupation weight. In an overlapping parent and daughter cell, the possible daughter momenta are the intersection of its cell with [p−q,p+q]; the other triangle inequality is automatic on these retained positive hard momenta. Bracket the two energy inequalities after splitting at the daughter-cell clipping transitions. Report the union in parent momentum separately from the sum of individual pair intervals, which can overlap near a piecewise cell boundary.

For the native prepared state, the smallest cold node q=1.1901e−6 has a parent range about [1.426752,1.430654], of width .00390204. The largest node q=8.9247e−6 covers about [1.414146,1.443407], of union width .02926085. For all eight modes in both hard reconstructions, initially and at the last recorded time, participating parent and daughter cells belong only to the two neighboring cells spanning [1.413889,1.466667]. Increasing the root-bracketing mesh from 16 to 64 preserves these cells and the union widths within 6e−13 absolute. The two cells independently emerge in the [occupation-substitution audit](INDEPENDENT-EARLY-POLE-INITIAL-COLD-SOURCE-HARD-CELLS-2026-11-27.md).

This verifies support in the finite dressed packet geometry, not the absence of contributions from altered or broadened dispersions, other hard supports, or outside the original cold scalar window. It does not assert that dynamics elsewhere on the hard grid are irrelevant.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
for n in 16 64; do
  python -W error early_pole_first_window_kinematic_support_reference.py \
    results/early-pole-first-window-64-memory-2026-11-27.json.gz --samples "$n"
done
```

Raw: `results/early-pole-first-window-kinematic-aperture-2026-11-27.json`.
