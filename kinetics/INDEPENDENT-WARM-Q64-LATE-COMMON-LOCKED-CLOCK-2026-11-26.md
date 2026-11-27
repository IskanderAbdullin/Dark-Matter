# Lock the stiff substeps while changing the late common pair chart order

The [late common-chart replay](INDEPENDENT-WARM-Q64-LATE-WHOLE-PAIR-REPLAY-2026-11-26.md)
could not pull apart angular/radial sensitivity and the
fallback in the inner implicit integration. Repeat both
orders on the exact same q64 parent, physical mesh, FRW
interval and old exterior gravity with **two mandated
equal prompt/elastic microsteps** of 0.001985628746 each
in Hubble units. Do not let the solver silently split them.
The stored common-chart end shapes, passively mapped back
to the midpoint, only initialize Newton.

Both orders accept both prescribed roots without switching
the initial-guess strategy, with scaled residuals below
8e-12 and the energy/number ledgers below 2.6e-14. The pair
kicks still use the default (3,7,6,2) or finer
(4,10,10,4) complete common chart; all else is matched.

| Comparison at the endpoint | Largest shell L1 over baseline represented moment |
| --- | ---: |
| Earlier adaptive solves, order changed | 9.859e-5 |
| Two equal inner steps, order changed | 2.439e-6 |
| Default order: locked versus earlier adaptive | 9.868e-5 |
| Finer order: locked versus earlier adaptive | 1.534e-14 |

The matched cold slopes are −5.13313820e-5 and −5.13315978e-5
per ln(a). The locked comparison therefore reduces the
earlier apparent full-shell order spread by about a factor
of 40. The finer locked result reproduces its adaptive
counterpart to ledger precision. This experiment separates
the two chosen numerical controls on **one** interval,
rather than asserting whole-trajectory convergence or
proving an error bound for either chart.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/warm-multishell-corrected-q64-stage-23-2026-11-26.json.gz
b=results/warm-multishell-corrected-q64-stage-24-2026-11-26.json.gz
c=results/warm-q64-late-full-common-default-2026-11-26.json.gz
d=results/warm-q64-late-full-common-fine-2026-11-26.json.gz
python -W error warm_q64_common_chart_locked_clock_reference.py EOS.dat "$a" "$b" "$c" locked-default.json
python -W error warm_q64_common_chart_locked_clock_reference.py EOS.dat "$a" "$b" "$d" locked-fine.json --fine
python -W error warm_q64_common_chart_locked_clock_audit_reference.py "$b" "$c" "$d" locked-default.json locked-fine.json
~~~

[Matched-clock audit](results/warm-q64-late-common-locked-compare-2026-11-26.json);
[default root and state](results/warm-q64-late-common-locked-default-2026-11-26.json.gz);
[finer root and state](results/warm-q64-late-common-locked-fine-2026-11-26.json.gz).
