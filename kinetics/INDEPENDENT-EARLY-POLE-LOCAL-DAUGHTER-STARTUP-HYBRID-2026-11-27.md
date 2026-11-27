# Refining just the first outer interval recovers 90% of the fixed-inner gap

Compare the same held native/local states over the original
first window at identical inner elastic width and 64 elastic
solves. Use four fine outer steps over the first sixteenth of
the exposure and fifteen original coarse outer steps after it.
This (19)-interval hybrid reproduces each saved cumulative
count in the first block exactly against the all-fine control.
Its individual elastic and original global ledgers pass.

| Outer schedule | Slow half-step count | Mean local minus native cold rate |
|---|---:|---:|
| Uniform 16 macro × 4 elastic | 32 | +2.12757801e-4 |
| Refine only first original interval, then coarse | 38 | +3.67384321e-4 |
| Uniform 64 macro × 1 elastic | 128 | +3.84570917e-4 |

The hybrid recovers 89.997% of the paired difference from
uniform coarse to uniform fine. Pasting the previously measured
starting-interval gap into the old coarse integral, without
evolving the hybrid, predicted +3.67381875e-4; the actual
hybrid differs by only +2.45e-9 in mean rate. This causally
tests the startup refinement within this held benchmark.
The residual versus uniform fine remains -1.71866e-5, and
the uniform fine result itself is not the unsplit reference.
Half-step counts describe work structure, not a wall-time
speedup or a general error bound.

The variable schedule is opt-in and uses fixed positive integer
weights specifying relative macro step lengths. The original
uniform default reproduces an independent held one-step
comparison byte-for-byte. Existing checkpoints include solver
identity; do not mix prior and edited solver checkpoints.

Use `early_pole_daughter_startup_hybrid.py` with the two
prepared cases and default `--startup-outer-intervals 1`.
Then run `early_pole_local_startup_hybrid_compare.py` on
the three native/local profile pairs.
