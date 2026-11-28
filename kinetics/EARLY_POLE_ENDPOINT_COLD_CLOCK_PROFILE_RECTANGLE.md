# All four choices of the final clock give almost the same cold profile gap

Use the exact 17 bins and four verified endpoints from the matched-clock profile receipt. On each normalized response, compare both same-clock S2/S4 pairs and both mixed-clock pairs. A triangle inequality bounds the total span of the four differences by the sum of the two intra-chart clock TVs. This bound concerns only the two calculations from each saved restart.

| Response | S2 own clock TV | S4 own clock TV | Span bound as fraction of lowest cross-chart TV |
| --- | ---: | ---: | ---: |
| Forward number | 8.25e-6 | 7.82e-6 | 0.0794% |
| Inverse number | 1.22e-5 | 1.19e-5 | 0.1192% |
| Forward kinetic excess | 1.14e-5 | 1.05e-5 | 0.1149% |

In all six matched comparisons (three channels at two clocks), the largest absolute signed common-bin contribution is still q=[3.30600e-6,3.97454e-6]. The four forward TVs range from 2.023732% to 2.024273%, even when the clocks are intentionally mismatched. No finer or earlier time-step error, absolute amplitude or physical convergence is bounded by this rectangle.

Use the same q-panel and ordered state arguments as `early_pole_endpoint_cold_matched_clock_profile.py`, but replace `--integrated-clock-receipt` with `--profile results/early-pole-cold-matched-clock-profile-2026-11-28.json` for `early_pole_endpoint_cold_clock_profile_rectangle.py`. The latter validates the exact four state hashes and both previously archived matched diagonals.
