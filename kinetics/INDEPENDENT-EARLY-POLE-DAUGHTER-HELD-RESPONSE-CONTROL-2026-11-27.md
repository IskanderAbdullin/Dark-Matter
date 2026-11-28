# Held unlike response does not remove the first-window daughter discrepancy

Use the same fine hard supports, exact prepared snapshots and cold
seed, duration 3.3333333333333335e-7, and 16 split macro steps.
Repeat the native and outside daughter preparations with the archived
fixed unlike dielectric rather than refreshing the response. Do this
with unlike-only and two-channel elastic operators. Also repeat the
aperture-local daughter preparation with both elastic channels. All
of the matched split variants pass the unchanged energy and entropy
gates.

| Daughter preparation and elastic channels | Active response, mean shift | Held response, mean shift |
|---|---:|---:|
| Outside, unlike only | -0.000752345519 | -0.000752342928 |
| Outside, both | -0.000741836280 | -0.000741834490 |
| Aperture local, both | -0.0000100299144 | -0.0000100295605 |

The two outside discrepancies from refreshing the response are below
2.6e-9 in mean source. In the aperture-local control the change from
holding the response is about 3.54e-10. The unsplit held-graph shift is
+0.00051724531, opposite to either 16-step split result, and its
cumulative effect crosses from negative to positive between saved
exposures 1.25e-7 and 1.4583e-7. The matched switched response therefore
does not account for the discrepancy at 16 steps. The remaining
difference is between the finite split and simultaneous time
integrations; this comparison alone does not distinguish the finite
width of the splitting from the finite implicit elastic step, or
give a bound at other widths or for active-screening dynamics in
general. A matched 64-step held control tests the next width.

From `kinetics/`, use the paired archived first-window bundle and the
prepared wrapper. For the outside and aperture-local runs select
`--case outside-daughter` or `--case local-daughter`; compare to
`--case native`. Use `--hold-screening` for the frozen response and
`--only-unlike` for the unlike-only operator. The two comparison
scripts check matched initial data and emit the recorded rates and
histories. No cosmological production or physical dominance follows.
