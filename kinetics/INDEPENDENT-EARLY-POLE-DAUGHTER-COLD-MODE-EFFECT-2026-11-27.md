# Cold-bin cancellation in the 16-step aperture-local daughter control

Resolve the cumulative cold-count shift on the eight existing cold bins
at each saved snapshot of the held-response 16-step run. Divide each
count difference by elapsed exposure. This uses the same native baseline
for the two independent daughter substitutions, with no new evolution
or change of cold mesh. The terminal sum agrees with the total scalar
ledger within the specified summation tolerance.

| Daughter control | Final signed mean shift | Sum of absolute mode shifts |
|---|---:|---:|
| Initially active two cells | -1.00296e-5 | 1.20903e-4 |
| Other 18 retained cells | -7.41834e-4 | 7.41834e-4 |

In the local control all eight mode shifts are negative in the first
three saved positive-time snapshots. At the end the two largest
opposing components are the fifth and sixth cold bins, about -4.83e-5
and +5.07e-5, respectively. Summing the eight modes cancels 91.7% of the
total absolute shift. For the outside preparation all mode shifts
remain negative at the saved times. The small local net shift in this
split control therefore does not mean that each mode is insensitive.
This calculation cannot establish the unsplit spectrum or cold mesh
convergence. Use `early_pole_daughter_cold_mode_effect_reference.py`
with the corresponding native and prepared frozen-screening outputs.

## Width refinement at fixed response

The matched [64-step held control](INDEPENDENT-EARLY-POLE-LOCAL-DAUGHTER-HELD-TIME64-2026-11-27.md)
has a positive terminal difference in all eight saved cold bins. The
fifth and sixth bins now give +1.51476e-4 and +1.62075e-4, while
the terminal sum is +3.84571e-4. Terminal cancellation as defined
here is zero, rather than 91.7%, and the terminal sum of magnitudes
increases by about 3.18 times. This is a change in the signed
distribution across the same cold bins, not merely a change in a
numerically subtracted scalar total. It does not show that all cold
modes have converged, or establish the spectrum of the simultaneous
integrator. The fifth and sixth bins are reported in one-based order.

The 64-step receipt uses the same analysis command with the accepted
held native and local prepared raw outputs at 64 steps.
