# A small count-clock drift does not imply small local collision-clock drift

Recover the fine D8 1024-step and third D8 1024-step endpoints from each chart's own archived .04 state. In both cases the recomputed coarse trace matches the independently stored charge/continuity/count records exactly. The coarse scalar gap shifts by about −1.09e−5 relative and the fitted photon energy by about −1e−11. Compute actual cold coefficients at k=.0015 on the two elapsed-time samplings, using the full separately refitted occupations, gap, photon energy and residue. Fine D8 is tested at q=3.5e−6; third D8 at q=3.5e−6 and 1.5e−5.

| Grid | q | Forward coarse/fine −1 | Inverse coarse/fine −1 | Forward shift divided by count drift, magnitude |
| --- | ---: | ---: | ---: | ---: |
| Fine D8 | 3.5e−6 | −0.40962% | −0.45041% | 125.61 |
| Third D8 | 3.5e−6 | −0.43294% | −0.48448% | 132.47 |
| Third D8 | 1.5e−5 | −0.31174% | −0.33713% | 95.38 |

The scalar counts move by about −3.26e−5 on either chart. The angular mesh error (<2.77e−13 relative) is distinct from the effect of changing the evolving step count. The fine-time rate for each benchmark exactly matches the independently published [fixed-k values](INDEPENDENT-EARLY-POLE-ENDPOINT-D8-HELD-K-RATE-POINTS-2026-11-28.md). Local cold rates can be over a hundred times more sensitive to the final time-step choice than the terminal count. These are resolved point-diagram comparisons, not a bound on how a whole band or whole history converges.

Reproduce the coarse states using `early_pole_endpoint_charge_cold_clock.py` and the relevant fine or third archived D8 .04 state and paired .05 trace, as in the [third-only receipt](INDEPENDENT-EARLY-POLE-ENDPOINT-THIRD-D8-COLD-CLOCK-2026-11-28.md). Then run `early_pole_endpoint_charge_cold_clock_points.py COARSE FINE OPTICAL_RECEIPT`, with `--q 3.5e-6` for the fine card and the default two-point set for third. `early_pole_endpoint_peer_d8_clock_points_receipt.py` checks the paired receipt SHA-256 values, coarse ratios and independent late fixed-k baselines.
