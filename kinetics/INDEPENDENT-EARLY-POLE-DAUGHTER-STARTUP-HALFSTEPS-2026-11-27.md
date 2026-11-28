# First-interval replay puts the discrepancy in the slow half-step sampling

Repeat only the first common interval (2.0833333333333335e-8)
of the held 16×4 and 64×1 native/local split profiles. Retain the
original full-duration initialization, fixed support, inner
elastic step, and energy budget. Capture the cold counts just
before and after each slow half-step. All four replays reproduce
the accepted original first-interval gains exactly and pass
the individual elastic energy gate. The recorded cold counts
are unchanged in between the slow solves, including across
elastic steps.

| Schedule over first common interval | Paired integrated local-minus-native cold gain |
|---|---:|
| 16×4, two slow half-steps | -1.77275718e-10 |
| 64×1, eight slow half-steps | -1.25734360e-10 |
| Difference | -5.15413580e-11 |

The first coarse half-step alone is -1.61567874e-10
and the second is -1.57078428e-11. The fine half-step
contributions begin at -4.0392e-11 and fall in magnitude
to -3.9269e-12. Thus the placement of slow updates changes
how much of the initial suppression is integrated before
the hard collision state relaxes. The elastic updates themselves
leave the tracked cold counts unchanged in this observed split.
This is an audit of the numerical stages, not an additive
attribution of a physical production rate.

From `kinetics/` run `early_pole_daughter_startup_halfstep_audit.py`
with the fine raw bundle and the original accepted native/local
pairs under `--outer16` and `--outer64`. The partial solver keeps
the archive's original full width and records two or eight
slow half-steps through the shared first interval only.
