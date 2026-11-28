# Thin charged faces sharpen the decay approximation but introduce fast scales

The [edge-padding test](INDEPENDENT-EARLY-POLE-EDGE-PADDED-CHARGED-PROJECTION-2026-11-28.md) recovers the direct frozen histogram integral by restricting interpolation to thin face cells. The probability factors also make the source nonlinear in those cells. Differentiate that source while holding the chart (masses, apertures, centroid weights and expansion factor) fixed.

For every event, differentiate its forward and inverse factors with respect to the scalar cell, two daughter cells and two parent cells. Weight the result by the same five signed deposition coefficients and divide by the destination volumes. This gives a sparse Jacobian, with left null vectors corresponding to the paired numbers and the fixed energy. The analytic directional derivative is independently checked with centered occupation increments. This is not the full Jacobian of a live scalar mass or evolving expansion.

| Archive | Unpadded inverse fastest charged diagonal | Face fraction 1/4096 |
| --- | ---: | ---: |
| Fine | 4.67e−2 | 2.83e−8 |
| Third | 7.61e−3 | 1.78e−8 |

The entries are measured in log-scale-factor units; inverse diagonal magnitudes are indicators of fast local response, not certified eigenvalue scales or stability bounds. The fastest charged response at the thinnest tested pad is in a daughter cell on either archive. The calculation warns against carrying the static 4e−5 error gate over to a time integration without explicit occupation and conservation controls. No time integrator is certified here.

Run `OPENBLAS_NUM_THREADS=1 python -W error early_pole_padded_stiffness_pilot.py` from `kinetics/` and compare the [receipt](results/early-pole-padded-stiffness-2026-11-28.json).
