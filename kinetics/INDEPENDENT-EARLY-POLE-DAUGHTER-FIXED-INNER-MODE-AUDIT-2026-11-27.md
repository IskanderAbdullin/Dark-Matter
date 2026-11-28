# All eight cold modes respond positively at the fixed refined elastic width

Resolve the [three held split profiles at the same elastic step](INDEPENDENT-EARLY-POLE-LOCAL-DAUGHTER-FIXED-INNER-MIDDLE-2026-11-27.md)
mode by mode. Each compares identical native and two-cell local
daughter states across profiles, with the same cold grid and seed.
Use only exact shared saved times and the same first-window
duration. Each profile passes the original global closure gates
and the explicit energy gate on every individual elastic step.

At the end of the window all eight local-minus-native mode
differences are positive in each of the 16×4, 32×2 and 64×1
profiles. Increasing the frequency of slow half-step interleaving
therefore changes the size, not the terminal sign, in this
matched-inner comparison. For the coarsest-to-finest comparison
the sum of absolute terminal mode errors is 1.71813e-4; the fourth
and fifth bins in one-based order account for 73.84% of it. The
intermediate-to-finest error is 3.79596e-5; the same bins account
for 74.61%. Each intermediate terminal mode is between the coarse
and fine terminal modes. These are rate shifts averaged over the
first window, not instantaneous or physical spectral estimates.

This does not assume mode convergence. It differs from the
[16×1 coarse result](INDEPENDENT-EARLY-POLE-DAUGHTER-COLD-MODE-EFFECT-2026-11-27.md),
where oppositely signed late mode effects nearly cancel in the
total. The 16×4 and 16×1 profiles use the same 16 outer intervals,
but the former uses a smaller elastic step.

The receipt is reproduced by
`early_pole_daughter_fixed_inner_mode_compare.py` using accepted
native/local raw pairs and the validated substep and mode-audit
functions.
