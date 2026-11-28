# Check heavy-parent recycling on the matched warm states

The [closed-count comparison](INDEPENDENT-PORTAL-CLOSED-NUMBER-BUDGET-2026-11-26.md)
requires an actual new source of count per visible entropy for
a dominant scalar component. Repeat the earlier
[heavy-parent rate diagnostic](INDEPENDENT-ARCHIVED-HEAVY-RECYCLING-2026-11-26.md)
on the actual q16 and q32 16-block endpoints at 10.442 MeV.
This is the same frozen optimistic formation rate, not a
time integration. Include both opposite and same charges
in the stated bound, omit final heavy Pauli blocking, and
credit every formed heavy particle with one net count after
its possible prompt decay. No inverse recycling is included.

| At 10.442 MeV | Scalar initial-state gross increase / (N H) | Stable opposite-charge increase / (N H) | Optimistic bound including both same-charge modes |
| --- | ---: | ---: | ---: |
| Evolved q16 | 3.060×10⁻⁷ | 5.725×10⁻⁶ | 1.748×10⁻⁵ |
| Evolved q32 | 3.881×10⁻⁷ | 6.747×10⁻⁶ | 2.063×10⁻⁵ |
| Same q32 state cropped to q16 | 3.089×10⁻⁷ | 5.792×10⁻⁶ | 1.768×10⁻⁵ |

The widened q32 state admits more endothermic reactions;
the same-state crop separates the instantaneous momentum
cut from the small difference in evolved inner shapes.
Free redshift on the retained stable support then closes
the geometric threshold at relative factors 1.550 and 3.100.
That statement does not apply with scattering or restored
older exterior populations. The rates use represented counts
and the visible radiation clock as in the earlier control
and are **not** an integrated number-multiplication bound.
They do not establish or exclude a separate source history,
or an unknown number-changing channel.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error archived_heavy_recycling_gate_reference.py EOS.dat \
  "$D"/warm-q{16,32}-10p442-heavy-input-shape-2026-11-26.json --crop-last 16
~~~

The reproducibility command uses descriptive archive filenames;
only the corresponding snapshot-name strings differ from the
compute record.
