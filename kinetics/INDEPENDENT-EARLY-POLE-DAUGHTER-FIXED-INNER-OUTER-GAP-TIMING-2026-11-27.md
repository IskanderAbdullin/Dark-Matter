# The fixed-inner outer-width cold-count gap accrues early

Compare the 16×4, 32×2 and 64×1 held split profiles at the
[same 64 elastic step widths](INDEPENDENT-EARLY-POLE-LOCAL-DAUGHTER-FIXED-INNER-MIDDLE-2026-11-27.md).
Run native and local swaps in each and subtract their cumulative
cold gains. The audit cross-checks the fixed physics, exact
initial hard/cold data, individual elastic receipts and snapshot
times. The outer 16 grid is the common count-reporting grid;
instantaneous comparisons use only times with saved hard states
in every run.

| Profile minus 64×1 | Terminal paired count gap | Fraction by first common interval | Fraction by first two |
|---|---:|---:|---:|
| 16×4 | -5.72710e-11 | 89.996% | 97.535% |
| 32×2 | -1.26532e-11 | 90.574% | 97.722% |

The first common interval is exposure 2.08333e-8, one sixteenth of
the window. By exposure 8.33333e-8 both profiles have accumulated
over 99.76% of their terminal discrepancy. Nevertheless the
instantaneous local-minus-native rates at the shared saved
states agree to 8.48e-9 and 1.65e-9, respectively. The largest
paired hard-moment component gaps on those same states are
2.77e-10 and 5.77e-11.

Thus the fixed-inner difference is mainly an early integrated
difference, not a persistent disagreement in the sampled rates.
This is consistent with a sensitivity to the placement of the
slow substeps during rapid startup. The saved endpoints do not
resolve transients inside each interval, so this does not
identify an exact mechanism or prove dense-output agreement.
The result is in solver exposure, not elapsed wall time.

Reproduce with `early_pole_daughter_outer_gap_timing_audit.py`
and the three accepted native/local raw pairs.
