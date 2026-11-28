# Measure how much nonmonotonic variation a positive frozen witness uses

A [per-cell monotone envelope](INDEPENDENT-EARLY-POLE-DAUGHTER-MONOTONE-ENVELOPE-2026-11-27.md)
gave only negative upper estimates. To relax monotonicity
without prescribing a pointwise slope, measure upward variation
(V_+(f)=sum_i(f_{i+1}-f_i)_+) for each continuous
piecewise-linear profile. Minimize the sum for both daughter
arms with the same exact modeled moments, first for feasible
moments alone and then for a positive paired rate margin.
This variation is a cumulative sum, not a fraction: it
can exceed one when a profile repeatedly rises.

| Hard depth × cold modes | Feasible witnesses, three moments | Positive witnesses, three moments | Positive witnesses, two moments |
|---|---:|---:|---:|
| 1 × 8 | .001897 | .142800 | .027479 |
| 1 × 32 | .001897 | .143456 | .027904 |
| 4 × 8 | .001996 | 1.520480 | .414789 |
| 4 × 32 | .001996 | 2.021108 | .407714 |

These are verified finite witnesses at 128 subdivisions.
For example, the depth-four/cold-32 value at 64 subdivisions
is 2.208696 instead. Matching both outer aperture values
raises the 128-subdivision value to 2.022351. The table is
not a claim of convergence to the minimum over continuous
physical distributions. The checked linear-program dual and
primal witnesses agree closely at each resolved finite mesh.

There is also a conservative numerical lower estimate that
does not rely on representing the occupation by those nodes.
The future supremum and past infimum of one occupation are both
nonincreasing. Each differs from it by at most its upward
variation on that cell. Apply either as a one-sided comparator
to the monotone cutoff dual residual; use the smaller integral
of its positive or negative part as the variation penalty.
Summing over both arms and cells, the largest penalty times
the total upward variation bounds the increase of the
paired objective over its monotone upper value. With exact
integrals this is a continuum inequality. The quadrature
evaluations give necessary lower estimates of .13530 (two
moments) and .13757 (three moments) for a nonnegative
depth-four/cold-32 contrast. There remains a substantial gap to
the explicit 128-subdivision positive witnesses. Neither a
cutoff prior nor arbitrary oscillation is inferred to be physical.

Use `early_pole_daughter_minimum_upward_variation.py` with
`--subdivision 64` and 128 for the full, matched-exterior,
and `--two-moments-only` controls. Run the monotone script at
12 and 48 cutoff subdivisions for the certificate and
combine with `early_pole_daughter_minimum_upward_variation_compare.py`.
The analytic two-node variation test is independent of the archive.
