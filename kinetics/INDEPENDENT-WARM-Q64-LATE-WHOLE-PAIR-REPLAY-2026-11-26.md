# A complete reversible pair chart on the final coarse q64 interval

The [earlier one-block common-chart replay](INDEPENDENT-WARM-PAIR-COMMON-CHART-2026-11-26.md)
was near the native cold spike. Repeat on the later 24th
q64 interval from 7.423463 to 7.394049 MeV. Replace the
complete pair operator in both directions. Keep the same
actual parent, old-exterior gravity, expansion clock and
selected prompt/elastic operators. Offer the accepted
nine-cell endpoint only as the Newton initial guess.
No q32 profile is imposed on the new collision targets.

| Accepted solution | Cold-cell change per ln(a) |
| --- | ---: |
| Nine-cell baseline | −5.07827452e-5 |
| Whole common chart, default (3,7,6,2) | −5.07826582e-5 |
| Whole common chart, finer (4,10,10,4) | −5.13315978e-5 |

The cold decline survives the switch to either full chart.
The largest whole-shell difference in number, energy and
pressure between the two full-chart outputs is 9.86e-5 of
the corresponding nine-cell baseline moment. The relative
cold-slope difference is about 1.1%. Both guided solves needed
the existing stiff-solver fallback; their last unsuccessful
guided residuals were 8.27e-12 and 1.07e-11. Their separate
collision/FRW energy and number ledgers close within 3.5e-14.
Since the fallback may have accepted different internal time
substeps, this comparison **does not isolate angular quadrature**
from the difference in stiff time updates. It does not certify
whole-domain convergence.

The old exterior records survive rebasing unchanged. All three
runs write four new positive receipts; their added light counts
are about 6.0–6.9e-28 of the retained light count. This is a
single selected interval with fixed archived initial data,
not a full alternative trajectory or a real unknown exterior.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/warm-multishell-corrected-q64-stage-23-2026-11-26.json.gz
b=results/warm-multishell-corrected-q64-stage-24-2026-11-26.json.gz
c=results/warm-q64-late-full-common-default-2026-11-26.json.gz
d=results/warm-q64-late-full-common-fine-2026-11-26.json.gz
python -W error warm_q64_late_common_pair_chart_audit_reference.py "$a" "$b" "$c" "$d"
~~~

For exact rerun of the two nonlinear solves, in a detached
working copy at fork commit `5ec7d34db4d5544be62e3c3b04b2baf36b92b6d8`,
copy the `warm_q64_late_common_pair_chart_reference.py` replay script
into its `kinetics` directory (preserving the pinned source kernel).
Invoke there with `EOS.dat "$a" "$b" default.json` and separately
with `EOS.dat "$a" "$b" fine.json --fine`. The prediagnostic
whole-chart source SHA-256 begins `721024d5361aae80`.

[Default accepted state](results/warm-q64-late-full-common-default-2026-11-26.json.gz);
[finer accepted state](results/warm-q64-late-full-common-fine-2026-11-26.json.gz);
[raw-state and outside ledger audit](results/warm-q64-late-full-common-compare-2026-11-26.json).
