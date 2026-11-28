# Combine nonoverlapping charged-pion slices and a separate neutral input

The earlier charged-pion results could not simply be summed:
the `|Fπ|²` threshold hypothesis overlapped the new
[BABAR moment bound](INDEPENDENT-WARM-BABAR-DISPERSIVE-GAP-2026-11-26.md).
This comparison first separates the charged spectrum into
three disjoint intervals, using:

- Below .30 GeV, either an explicit `|Fπ|²` cap on the
  vacuum two-pion phase or the [positive theory moment](INDEPENDENT-WARM-DISPERSIVE-THEORY-PION-THRESHOLD-2026-11-26.md).
  Never apply both in the same case.
- From .30 to √.10 GeV, the larger new BABAR allowance
  on that interval only. Assign none of it to KLOE.
- From √.10 through the selected upper edge, the existing
  KLOE measured-band model and bounded massive-pion Bose
  factors. No KLOE/BABAR fit or statistical combination.

The neutral `π⁰γ` [conditional moment](INDEPENDENT-WARM-NEUTRAL-PION-PHOTON-BUDGET-2026-11-26.md)
is a different channel. Use it only below .60 GeV with
the stated timelike cap on `Fhat`. Charged pion heat uses
its separate spacelike cap. Calculate the matched charged
baseline at its own cooling history; credit the BABAR,
neutral, and (where used) theoretical threshold allowances
in full at the start. Evaluate the actual positive quadratic
upper after those deposits and possible direct decay count,
rather than using the looser secant approximation. No
other added hadronic or scalar source is assumed.

| 50 MeV vacuum comparison | Finer scalar Ωh² upper |
| --- | ---: |
| Bose-weighted KLOE with only the additional BABAR slice | .0306830 |
| Plus sub-.30 timelike cap 1, elastic cap 1, and neutral amplitude cap 1 GeV⁻¹ | .0308158 |
| Plus sub-.30 timelike cap 4, elastic cap 4, and neutral amplitude cap 10 GeV⁻¹ | .0356479 |
| Replace the timelike threshold cap with the dispersive theory allowance, keep wider elastic and neutral caps | .0418371 |

No result above predicts a late relic or uses data to establish
the proposed form-factor caps, thermal response, or modes not
included here. The last, more conservative row still has
5.794×10⁻¹⁰ GeV⁴ of positive comoving input headroom before its
upper could reach .060. A change in the number per energy
of extra particles (in particular external light scalar
injection) cannot use this same headroom.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_joint_babar_pion_budget_reference.py EOS.dat KLOExSec.data \
 "$D"/warm-energy-weighted-recycling-kernels-2026-11-26.json \
 "$D"/warm-babar-dispersive-gap-2026-11-26.json \
 "$D"/warm-neutral-pion-photon-budget-2026-11-26.json \
 "$D"/warm-dispersive-theory-pion-threshold-2026-11-26.json
~~~
