# Adaptive stepping for the bounded local charge-face replay

Use the same archived start, face fraction 1/4096 and 5e−5 interval from the [two-width fixed-time comparison](INDEPENDENT-EARLY-POLE-BOUNDED-FACE-LONG-WIDTH-2026-11-28.md). Compare one bounded midpoint step with two half steps. The controller checks separate volume-weighted L1 occupation differences on the three species, plus a pointwise charged occupation difference. Accept only the two half-step state; do not extrapolate beyond the admissible occupations. Rejected candidates never contribute to the history or its energy budget.

| Archive | Accept/reject, error scale .002 | Accept/reject, error scale .0005 | Relative difference to 512 fixed midpoint steps (strict run) |
| --- | ---: | ---: | ---: |
| Fine | 45 / 2 | 91 / 9 | −5.27e−6 |
| Third | 71 / 12 | 132 / 20 | −5.15e−6 |

The more permissive runs agree with the fixed reference to 2.18e−5 or better. Each attempted macro step evaluates three candidate midpoint advances, so counts alone overstate the saving relative to a fixed midpoint step. The recorded bounds hold; the independent energy and paired number residuals are retained. This is a time integrator comparison at only one face width. It does not remove the much larger width dependence or certify the full history.

Run `OPENBLAS_NUM_THREADS=1 python -W error early_pole_adaptive_face_extents_pilot.py GRID` with fine and third. The [combined receipt](results/early-pole-adaptive-face-extents-2026-11-28.json) adds the relative differences to the previous fixed midpoint results.
