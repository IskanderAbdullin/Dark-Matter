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

## Which static choices give curvature?

Replay the full packet cold-source computation exactly at zero,
half, and full replacement. Hold the initial parent occupation and
scalar seed as before. In one diagnostic vary the reconstructed
daughter occupations while freezing packet geometry at native. In
another move the packet geometry while freezing native occupation
logits. These diagnostics intentionally separate coupled model
inputs; neither is a new physical source or time-evolution kernel.

| Hard depth × cold modes | Full midpoint curvature | Reconstructed occupations at fixed geometry | Geometry only, fixed logits |
|---|---:|---:|---:|
| 1 × 8 | +3.29827e-4 | +3.29862e-4 | -6.9234e-9 |
| 1 × 32 | +3.46261e-4 | +3.46215e-4 | -8.5078e-9 |
| 4 × 8 | +1.27911e-6 | +1.19711e-6 | -1.0006e-10 |
| 4 × 32 | +1.61783e-6 | +1.46489e-6 | -1.6215e-10 |

The occupation reconstruction alone accounts for nearly all
the coarse curvature in this diagnostic. On the refined hard
mesh its contribution is smaller in absolute terms and is about
91–94% of the residual total. The difference between the full
curve and the two individual diagnostics is a mixed term; the
moving geometry with native logits alone is small. At fixed
geometry, replacing reconstructed midpoint occupations by the
pointwise average of their endpoints gives exact affine
cold-source interpolation up to 6.94e-18. This follows because
the gain-minus-loss packet expression is affine in the daughter
occupation with parent and scalar occupations fixed. It does
not imply that the reconstructed midpoint moment state agrees
with that pointwise average.

Use `early_pole_daughter_local_fraction_curvature_attribution.py`
on the same four fixtures. Raw per-mode curvature and mixed
terms are in `results/early-pole-daughter-local-fraction-curvature-2026-11-27.json`.
