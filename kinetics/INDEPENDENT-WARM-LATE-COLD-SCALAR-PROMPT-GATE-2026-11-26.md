# Test the prompt supply of genuinely low-momentum scalars at 7.394 MeV

The high-energy excess in the near-cut tests is not the
population of scalars with p/T below .1. At the corrected
q32 and q64 late archives, integrate the scalar daughter
count and the separately signed prompt gain/loss at six
physical p/T cuts. The parent is stationary with daughters
held fixed, and its allowed energy is bounded by the
sum of maximum represented incoming daughter energies.
Distinguish closed decay from decay including outgoing
daughters beyond the archive; both use a vacuum incoming
exterior in the folded diagnostic.

At p/T=.1 the initial scalar fractions are approximately
1.8755e−4 of the represented scalar number in both archives.
The folded direct prompt gain and loss there are both
zero in all independent orders. The exact two-body
minimum parent energy for those scalar momenta is
3.24662 GeV, above the incoming parent energy caps
of .43091 and .86043 GeV. A cold scalar with p/T≤.1
also needs a light inverse partner above 3.24582 GeV
in momentum. Even the separate 1.05 GeV pad of the
near-cut shape test cannot access this channel **directly**
with represented incoming daughters.

At p/T=1, where the minimum parent energy drops to
.33229 GeV, outgoing phase is important in q32:
the closed gain/loss are zero, but the open gain is
5.78e−12 times the represented scalar number per Hubble
time. A closed event for that cut requires a light momentum
at least .32428 GeV, above q32's .21499 GeV pad. Q64 can
access it and the closed net is +7.07e−13 in the same
units. At p/T=3 the net is about +6.085 in both archived
phases. The integrated scalar birth and loss still agree,
since this operator only redistributes particle momenta.

The geometric zero is conditional on the incoming support;
the archives have unresolved positive exterior cohorts,
which are not estimated here. Nor does this rule out
population changes from other reactions, expansion,
earlier source histories, or a preexisting heavy population.
The numbers at accessible cuts are local frozen currents,
not a proposal to integrate a rate for one Hubble time.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_late_cold_scalar_prompt_partition_reference.py EOS.dat \
 results/warm-corrected-q32-time48-step48-2026-11-26.json.gz \
 results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
~~~

[Angular and number controls, including exact recoil gates](results/warm-late-cold-scalar-prompt-partition-2026-11-26.json).
