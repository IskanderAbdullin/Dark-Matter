# The aperture-local daughter sign switches when the held split is refined

At the original hard supports, held elastic weights and common
3.3333333333333335e-7 exposure, evolve identical native and aperture-local
prepared inputs with 64 rather than 16 outer split steps. At each width
the initial hard and cold data match their corresponding archived
controls exactly. Both 64-step calculations pass the original
conservation requirements, including the local 5e-15 energy gate.

| Time method | Mean shift (local daughter minus native) |
|---|---:|
| Backward elastic split, 16 steps, active response | -1.00299144e-5 |
| Backward elastic split, 16 steps, held response | -1.00295605e-5 |
| Backward elastic split, 64 steps, held response | +3.84570917e-4 |
| Simultaneous charged-logit solve, held response | +5.17245309e-4 |

At fixed response the sign changes with time resolution. The 64-step
error in the mean shift relative to the simultaneous calculation is
about one quarter of the 16-step error. The two-width linear
cancellation gives +5.16104409e-4, leaving 1.14e-6 in the mean shift
relative to the simultaneous result. That consistency does not prove
formal order or a general error bound. Across the shared saved time
grid the largest cumulative difference drops from 1.76e-10 to
4.44e-11, or to 5.26e-12 after the same extrapolation. The integrated
effect crosses from negative to positive later in the 64-step
path (between 1.67e-7 and 1.875e-7 on the shared saved grid) than in
the simultaneous history (between 1.25e-7 and 1.46e-7).

This establishes the numerical sensitivity of the 16-step sign in
the specified fixed-response model. It does not separate the
error from alternating collision flows from the error of the
backward elastic step. An elastic-only subcycle at the original
outer width can distinguish those controls.

From `kinetics/`, run the prepared swap wrapper with
`--case local-daughter --hold-screening --steps 64` on the paired
n216 archive, and verify against the matched native output using
`early_pole_local_daughter_time_control_reference.py --held64`.
The receipt includes all earlier and new time profiles.
