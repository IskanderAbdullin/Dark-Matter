# Where does the vacuum upper stop excluding half of .120?

Change only the instantaneous reheating temperature in the
[from-empty vacuum charged-lepton calculation](INDEPENDENT-WARM-FROM-EMPTY-RECYCLING-ENVELOPE-2026-11-26.md).
Keep the central tabulated visible entropy/radiation history,
the same energy-weighted kernel and the same dark channels.
Use the positive invariant and time panels twice with finer
resolution on the second comparison.

| Start, MeV | Coarser Ωφh² upper | Finer Ωφh² upper |
| --- | ---: | ---: |
| 54.5 | .051713 | .049365 |
| 55.0 | .056492 | .053898 |
| 55.5 | .061620 | .058799 |
| 56.0 | .067177 | .064096 |

The refined result still rules out half of .120 under the
specified assumptions at 55.5 MeV; a coarser enclosure there
is inconclusive. At 56.0 MeV neither enclosure establishes
exclusion. This is a **numerical decision bracket for this
upper comparison**, not a physical transition, a fit,
a demonstration of half above 55.5 MeV, or an estimate of
a unique transition temperature. No extrapolated envelope or
extra hadronic input is included. The
[separate 50 MeV selected pion-band test](INDEPENDENT-WARM-PION-BAND-AND-SOURCE-HEADROOM-2026-11-26.md)
must not be used as the source for these temperature points.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_from_empty_temperature_sensitivity_reference.py EOS.dat \
  "$D"/warm-energy-weighted-recycling-kernels-2026-11-26.json
~~~
