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
