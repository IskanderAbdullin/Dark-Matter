# A fixed elastic orbit predicts which startup blocks to refine

At fixed 64 elastic leaves, compare a four-leaf macro block
with four one-leaf macro blocks. In the tangent at zero slow
coupling, the original order uses the raw cold-source samples
at the two endpoints with weight 1/2 each. The four-piece
trapezoid weights the endpoints by 1/8 each and the three
interior samples by 1/4 each. The difference, integrated over
the block, is additive across the 16 original macro blocks.
See the [fixed-orbit construction](INDEPENDENT-EARLY-POLE-DAUGHTER-FIXED-ORBIT-QUADRATURE-2026-11-27.md).

The two initial blocks score highest whether ranked by the
native arm alone, the local-daughter arm alone, the paired
cold contrast, or the sum of absolute per-mode corrections
in the two arms. At a budget of one, two, three, or four fully
refined blocks, enumerating choices minimizes the absolute
paired tangent residual by choosing the corresponding first
blocks. Later per-block corrections change sign; adding
refinement blocks need not reduce this particular signed
observable monotonically. No claim is made for another orbit.

| Placement | Tangent predicts fraction of paired gap recovered | Previously measured full slow-coupled fraction |
|---|---:|---:|
| First block only | 89.9968% | 89.9969% |
| First two blocks | 97.5359% | 97.5359% |
| Second block only | 7.5392% | not used for this prediction |

Adding the tangent corrections to the full coarse mean gives
predictions within 2.4e-9 of both accepted startup hybrid
means. This is an offline calculation from all fixed-orbit
samples; it does not establish a cheap online estimator.

## Predeclared placement control

At the same 19 outer intervals and 64 elastic leaves as the
first-only hybrid, refine **only the second** original block.
The coarse full mean plus the exact tangent correction predicts
a paired mean of +2.25730831e-4. Audit the complete accepted
histories and original energy and entropy gates. Before opening
the new arm receipts, require the placement correction to be
positive, the recovered full gap to be between 6.5% and 8.5%,
and the absolute mean prediction error to be at most 1e-7.
If those fail, report the discrepancy rather than treating the
offline tangent as a calibrated proxy.

From `kinetics/`, run
`early_pole_daughter_fixed_orbit_startup_allocation.py` with
the fixed-orbit, dyadic-order, one-block hybrid, and two-block
hybrid JSONs in that order. The higher momentum resolution,
response refresh, and physical relic inference are outside
this first-window method test.

## Result of the placement control

The second-only hybrid gives +2.257107313552442e-4 as the
paired final mean. It recovers 7.5390% of the full coarse-to-fine
difference, versus 89.9969% for the first-only hybrid. The
difference from the predeclared second-only prediction is
-1.50104e-10, and all three predeclared conditions pass. Both
placements use 19 outer intervals and 64 fixed-width elastic
solves; slow RHS calls are 616 and 568 per arm, respectively.
The second-only first coarse sample matches the full original
coarse path exactly.

Reproduce with `early_pole_daughter_startup_hybrid.py` using
`--refined-blocks 2` for native and local-daughter on the same
archive. Gate the results using
`early_pole_daughter_startup_placement_compare.py`. This confirms
placement sensitivity for this window, not universal timing gains
or a wall-clock comparison.
