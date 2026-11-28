# How much heavy recycling remains on the unmodified warm spectra?

Start from the [matched 10.442 MeV q16 and q32 snapshots](INDEPENDENT-WARM-10P442-HEAVY-RECYCLING-2026-11-26.md).
For an expansion factor `b` relative to that anchor, redshift each
represented occupation as `f_b(p) = f_1(b p)` and shrink the same
shell edges by `1/b`. This preserves comoving number without
fitting another intra-shell shape. Integrate the two positive
formation rates, multiplied by `2 b³/(N₁ Hrad)`, in `d ln b`
from the anchor until the retained support crosses the
respective pair thresholds. Here the 2 counts both newly
formed parents, each allowed to decay promptly into two dark
daughters and thus add one net dark particle. Use the conditioned
visible radiation clock, neglect final Pauli blocking, and
allow the conservative same-charge contribution used in the
fixed-state reference. The denominator `N₁` is the represented
total at the anchor, held at the **full q32** value for its crop.

| Retained input | Scalar-initiated gain / N₁ | Opposite-charge stable-initiated gain / N₁ | Optimistic total including same charges | Share before b=1.12 |
| --- | ---: | ---: | ---: | ---: |
| q16 | 1.3513×10⁻⁸ | 2.7785×10⁻⁷ | 8.4707×10⁻⁷ | 92.26% |
| q32 | 1.8903×10⁻⁸ | 3.5359×10⁻⁷ | 1.0797×10⁻⁶ | 89.96% |
| Same q32 state cropped to q16 | 1.3648×10⁻⁸ | 2.8127×10⁻⁷ | 8.5746×10⁻⁷ | 92.25% |

The added support raises this passive exposure by 27.46% between
evolved q16 and q32. Cropping the *same* q32 state removes
3.23×10⁻⁵ of the represented number but 20.58% of the calculated
exposure. More than 89% of each exposure accrues by b=1.12,
well before the last retained mode redshifts below threshold.

The last two time/phase-space quadratures (one no longer splits
outgoing angles at shell edges) agree to 0.077% or better in the
optimistic totals. The b=1 rates agree with the separate
fixed-state fold to 0.069% or better. The comoving number and
second momentum moment of the discrete passive remap are
preserved to 4.5×10⁻¹⁶ in this independent quadrature check.

**Scope:** this gives about one optimistically credited particle
per million initially represented, *on these passive spectra*.
It adds time information to the instantaneous rate test, not a
general limit on number multiplication. Positive rates were
sampled on the unmodified occupations; neither outgoing products,
elastic repopulation, inverse reactions, older unrecorded tails,
nor additional formation channels enter the integral. A photon
clock with additional positive densities would be faster at
given `b` for the same visible entropy history, but this does not
justify treating the passive rates as an upper bound with feedback.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_ballistic_heavy_exposure_reference.py EOS.dat \
 "$D"/warm-q{16,32}-10p442-heavy-input-shape-2026-11-26.json \
 --crop-last 16 > exposure.json
python -W error warm_ballistic_heavy_exposure_verify_reference.py \
 exposure.json "$D"/warm-early-10p442-heavy-recycling-2026-11-26.json \
 "$D"/warm-q{16,32}-10p442-heavy-input-shape-2026-11-26.json
~~~
