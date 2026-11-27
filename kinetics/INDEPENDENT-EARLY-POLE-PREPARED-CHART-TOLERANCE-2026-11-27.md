# Tenfold tighter integration does not remove the held-chart daughter interaction

Repeat the full four-case [unsplit input comparison](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-UNSPLIT-DAUGHTER-ADDITIVITY-2026-11-27.md)
with relative tolerance 3e-10 instead of 3e-9 and base absolute
tolerance 2e-15 instead of 2e-14. Keep the original support,
cold seed, initial prepared states, held graph and all reporting
times fixed. Each of the eight trajectories passes the original
closure gates. The largest change across reported count histories
is 1.05e-17.

| Prepared shift from native | Nominal | Tighter |
|---|---:|---:|
| Aperture local | +5.17245309084e-4 | +5.17245319492e-4 |
| Outside aperture | -7.71660589594e-4 | -7.71660568777e-4 |
| Entire daughter | -2.54257750271e-4 | -2.54257739862e-4 |
| Entire minus local minus outside | +1.57530238853e-7 | +1.57509422172e-7 |

The signed interaction persists: its magnitude changes by only
0.0133%. The finite local and outside shifts are insensitive to
this particular numerical tolerance change. The test does not
bound hard support, quadrature, physical kernel or other sources
of error; the interaction remains a property of this manufactured
input comparison.

Reproduce the four tighter cases as previously, adding
`--rtol 3e-10 --base-atol 2e-15` to the chart command. Feed
the four nominal and four tighter raw receipts into
`early_pole_prepared_chart_tolerance_comparison.py` in native,
local, outside and entire order.
