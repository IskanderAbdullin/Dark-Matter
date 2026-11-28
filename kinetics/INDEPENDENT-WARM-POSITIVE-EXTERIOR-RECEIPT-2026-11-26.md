# Measure positive exterior creation, not repeated exterior reaction use

The [gross peripheral criterion](INDEPENDENT-WARM-HEAVY-PERIPHERAL-GROSS-BUDGET-2026-11-26.md)
needs energy **received from** outside daughters by events
creating high heavy particles. The current one-way exterior
solver instead stores positive energy **sent to** newly
projected exterior cohorts. These two ledgers must not be
identified. In this step measure the second quantity and
test only the explicitly conditional single-use case.

Use the [four matched saved blocks through b=1.02](INDEPENDENT-WARM-LIGHT-OCCUPATION-NEXT-MATCHED-STAGE-2026-11-26.md).
Check that each previously stored positive cohort keeps
the same injection boundary, species, initial mean
energy and comoving count when rebased. The new cohort
receipt is the injected energy per particle times the
new cohort number, expressed in original-anchor comoving
volume units, without later redshift. The independently
tracked radiation-style energy injection has an additional
factor of the cohort's injection expansion normalized
to the original anchor; use that identity as a second
check. Only count the upper-momentum exterior receipts
in the new sum. There are 16 positive new high cohorts
for each support over this short continuation.

| At b=1.02 after 10.442 MeV | q16 | q32 |
| --- | ---: | ---: |
| New one-way high exterior creation, GeV⁴ | 1.76973×10⁻¹⁷ | 1.74670×10⁻²³ |
| New creation / full available energy upper | 3.97621×10⁻⁶ | 3.92492×10⁻¹² |
| One-use pool including the older pressure-only energy upper, GeV⁴ | 3.18781×10⁻¹⁵ | 7.83835×10⁻²¹ |
| This conditional pool / permitted gross-use budget | .0443222 | 1.25149×10⁻⁷ |

The one-use pool includes the known initial tracked
exterior energy and the independent older energy upper
without restoring those old momenta. It could cap
reactive use **only** if each particle can spend its
energy once, without re-entry followed by regeneration
or additional heating, and only with the documented
source history. No such use is simulated here. A small
recorded one-way receipt does not imply small cumulative
positive input if recycling is restored. Nor is this
short receipt the future source over a long trajectory.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_positive_exterior_receipt_reference.py \
 "$D"/warm-signed-fermion-with-charged-pions-2026-11-26.json \
 "$D"/warm-heavy-peripheral-gross-budget-2026-11-26.json \
 "$D"/warm-known-exterior-recycling-envelope-2026-11-26.json \
 "$D"/warm-early-q{16,32}-time32-checkpoint-2026-11-26.json.gz \
 "$D"/warm-light-occupation-next-q{16,32}-tape-2026-11-26.json.gz
~~~
