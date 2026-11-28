# Five hard resolutions against direct archived daughter shapes

Refine the hard aperture by factors 1, 2, 4, 8, and 16,
holding the same selected physical daughter interval, parent
occupation, cold number grid, and model equations within each
frozen fixture. Regenerate the published 8-cold-mode depth-1
and depth-4 initial fields exactly before making new initial
fixtures. Only evaluate the initial paired source in the
native geometry of each fixture; no new time trajectory is
created by this calculation.

| Hard depth | Fitted minus direct, 8 cold modes | Fitted minus direct, 32 cold modes |
|---:|---:|---:|
| 1 | -3.85839e-3 | -4.03580e-3 |
| 2 | +2.42059e-4 | +2.36720e-4 |
| 4 | -1.82804e-5 | -2.07398e-5 |
| 8 | +9.52396e-7 | -8.37872e-8 |
| 16 | +4.92633e-8 | +1.11995e-7 |

The signed discrepancy does not converge monotonically. In
particular, at 32 cold modes the net cancellation is stronger
at depth 8 than at depth 16. The sum of absolute mode-cell
differences, which does not let those errors cancel, decreases
at every refinement for both cold supports: it is 7.63e-7 and
8.87e-7 at depth 16, down from 3.86e-3 and 4.04e-3. This is
still a finite mesh test rather than an asymptotic order claim.
The direct archive calculation
varies by at most 1.25e-7 across all five hard depths for fixed
cold support; this still retains the archived shape prescriptions
and native dressed geometry rather than determining a physical
continuum. The two cold supports continue to give distinct
initial rate estimates.

Run `early_pole_daughter_pointwise_hard_ladder.py ORIGINAL DEPTH4`
on the original 8-mode archive and exact depth-4 8-mode
initializer from the [same-interval test](INDEPENDENT-EARLY-POLE-DAUGHTER-ORIGINAL-INTERVAL-HARD-REFINEMENT-2026-11-27.md).
All values are in the JSON receipt.
