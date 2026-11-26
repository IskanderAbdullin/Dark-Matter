# Carry the measured one-way exterior ledger across a plan change

The [original receipt](INDEPENDENT-WARM-POSITIVE-EXTERIOR-RECEIPT-2026-11-26.md)
covered four blocks. Append the two [instrumented blocks](INDEPENDENT-WARM-OBSERVED-SUBOPERATOR-OCCUPATION-2026-11-26.md)
after verifying their actual source, plan and checkpoint
digests. The second stage was initialized with its own
local expansion origin. Rebase each cohort back to the
first 10.442 MeV anchor using its original comoving
number and injection scale. Earlier cohorts persist;
they are never counted as freshly injected by the
second driver.

At relative b=1.0302 there are 24 positive high
cohorts created since the 10.442 MeV anchor for each
support. The exact integrated high receipt is
(2.64627×10^{-17}) GeV⁴ on q16 and (2.60662×10^{-23})
GeV⁴ on q32, adding 49.53% and 49.23% to the preceding
two-percent interval receipts. These are
(5.94560×10^{-6}) and (5.85722×10^{-12}) of the
earlier gross energy upper, respectively.

With older unknown shapes represented only by their
earlier conditioned upper energy, the hypothetical
one-use pools are .044444 and (1.25286×10^{-7})
times the partial gross-reaction budget. This does not
measure any reacting use: every exterior cohort in
this numerical model stays collisionless by
construction. In particular the small pool cannot
be used as a bound on energy spent multiple times
by an exterior-aware reaction solver.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_exterior_cross_stage_receipt_reference.py \
 "$D"/warm-signed-fermion-with-charged-pions-2026-11-26.json \
 "$D"/warm-heavy-peripheral-gross-budget-2026-11-26.json \
 "$D"/warm-known-exterior-recycling-envelope-2026-11-26.json \
 --parents "$D"/warm-early-q{16,32}-time32-checkpoint-2026-11-26.json.gz \
 --earlier "$D"/warm-light-occupation-next-q{16,32}-tape-2026-11-26.json.gz \
 --later "$D"/warm-observed-suboperator-q{16,32}-tape-2026-11-26.json.gz
~~~
