# First-window finite sensitivity to elastic radial quadrature

At 16 macro steps over 3.3333333333333335e-7, double the radial
quadrature order of both elastic packet graphs from (1, 1, 2) to
(2, 1, 2), leaving the initial hard moments, cold seed, aperture,
screening rule and other packet orders unchanged in each preparation.
The mean cold source is split count gain divided by the window length.

| Preparation | Default | Radial 2 | Relative change |
|---|---:|---:|---:|
| Native | 0.06694243427 | 0.06700051953 | +867.69 ppm |
| Cell-conservative | 0.06528878623 | 0.06543588591 | +2253.06 ppm |

Continuous below native, normalized by native, changes from 2.47025%
to 2.33526%. Initial mode rates match exactly between quadratures, so
the difference is in their subsequent joined evolution. Largest
refined elastic step defect is 2.0e-15 with a 5.0e-15 per-step
budget. Terminal absolute energy drift is below 4.9e-15, and the
two entropy closure magnitudes are below 1.7e-12.

The receipt checks matched controls and accepted conservation using
`early_pole_first_window_finite_elastic_refinement_compare.py` and
the two default and two refined raw split outputs. Reproduce the
refined legs with `early_pole_joined_elastic_packet_refinement_reference.py`
and `--projection steps` or `--projection cell-conservative`.
Two radial orders do not establish convergence; both elastic
graphs changed together. This test leaves transfer, angle, time,
and hard-aperture refinement separate.
