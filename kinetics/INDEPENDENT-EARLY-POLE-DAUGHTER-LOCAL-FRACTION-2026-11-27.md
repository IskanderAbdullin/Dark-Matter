# The initial aperture response is closer to linear on the finer hard grid

Interpolate the daughter's two hard moments between the original
and alternate archive in just the exact local momentum interval.
The charged state outside it and the scalar seed stay fixed. A
convex moment interpolation is an admissible manufactured control;
it is not a dynamical path or a claim that the occupations interpolate
linearly. Reconstruct and dress the occupations at each interpolation.

| Hard aperture depth | Cold modes | Half replacement: initial rate difference | Midpoint departure relative to half the complete shift |
|---:|---:|---:|---:|
| 1 | 8 | -0.0074254540 | +4.2529% |
| 1 | 32 | -0.0077790386 | +4.2615% |
| 4 | 8 | -0.0058341420 | +0.02192% |
| 4 | 32 | -0.0061162981 | +0.02644% |

Here midpoint departure is (half replacement minus half the complete
shift) divided by the magnitude of half the complete shift. The
eight- and 32-mode results both show much smaller curvature when
the hard aperture is refined. Each retained mode decreases over the
nine sampled fractions in each fixture; this is a finite sampled
statement, not a proof over all possible fractions or occupations.
The full local time response [survives as a different diagnostic](INDEPENDENT-EARLY-POLE-DAUGHTER-ORIGINAL-INTERVAL-HARD-REFINEMENT-2026-11-27.md).
This static near-linearity does not determine any finite intermediate
or extrapolated response.

From `kinetics/`, run `early_pole_daughter_local_fraction_static.py`
on the original 8-mode archive, the exact depth-4 8-mode initializer,
and their 32-cold-mode initializers, in that order. Generate the
original-hard 32-mode initializer with
`early_pole_rebin_held_initial_cold_support.py INPUT OUTPUT`;
the depth-4 version is retained in `results/`. This tool reconstructs
only the initial cold seed and source on the nested 32-mode grid.
It verifies exact hard data, exact old initial sources and grouped
cold counts, and does not synthesize trajectory histories.
