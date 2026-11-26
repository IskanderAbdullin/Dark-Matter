# Compare .52 and .58 GeV moment-matched hard receipts in expansion

Repeat the [completed .58 GeV stress](INDEPENDENT-WARM-COUPLED-COHORT-STRESS-2026-11-26.md)
with the same 132 archived exterior cohorts and mean birth energies,
but put the distinguished hard packet at .52 GeV. Compare the two
genuinely evolved q64 runs and the existing vacuum-padded control
on the identical two-interval clock (10.391 to 10.380162 MeV).

| Observable | .52 GeV | .58 GeV |
| --- | ---: | ---: |
| Initial bare high inverse energy rate / retained energy / H | 1.4915e-20 | 6.6391e-21 |
| Terminal bare high inverse energy rate / retained energy / H | 3.7792e-32 | 3.7792e-32 |
| End / initial inverse rate | 2.534e-12 | 5.692e-12 |
| New scalar exterior count / terminal vacuum-control scalar count | 1.274e-22 | 1.608e-19 |
| New light exterior count / terminal vacuum-control light count | 1.927e-22 | 3.452e-19 |

Changing where the hard packet was placed changes the initial high
inverse rate by a factor 2.25 and the tiny new exterior counts in
the saved ledger by more than three orders of magnitude. Nonetheless
the two final **retained** spectra agree at the numerical closure
scale: number, full-energy and pressure shell L1 differences
normalized to the vacuum-padded moments are below 9.2e-15. Both
cases still differ from a vacuum initialization in the retained
counts, at order 1e-5 (scalar) and 4e-6 (light). The computed
end-to-end difference in the inverse rate is only 2.5e-11 relative;
this is not an error or accuracy guarantee.

The .52 GeV restart checks the original moment constraints, exact
terminal hash, matched moment mesh, inclusive count (7.11e-15)
and segment energy closure (2.08e-15 or better). This is a
two-realization sensitivity experiment on the selected coupled
operators. The order-of-magnitude ledger changes are measurements on
these two runs; time refinement would be needed
to claim quantitative accuracy for these especially small receipts.
Neither case bounds other realizations, untracked incoming states,
or a dynamically complete heavy spectrum.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p case052/q64
zcat results/warm-matched-explicit-support-q64-2-blocks-2026-11-26.json.gz >case052/vacuum.json
zcat results/warm-coupled-two-packet-q64-hard052-2026-11-26.json.gz >case052/q64/checkpoint.json
python -W error warm_coupled_two_packet_completed_reference.py EOS.dat case052 case052/vacuum.json --hard .52
python -W error warm_coupled_two_packet_hardness_comparison_reference.py \
  results/warm-coupled-two-packet-q64-hard052-2026-11-26.json.gz \
  results/warm-coupled-two-packet-q64-hard058-2026-11-26.json.gz case052/vacuum.json \
  results/warm-coupled-two-packet-hard052-2026-11-26.json \
  results/warm-coupled-two-packet-hard058-2026-11-26.json
# Fresh solve in a separate directory:
MALLOC_ARENA_MAX=2 python -W error warm_cohort_two_packet_stress_reference.py EOS.dat rerun052 case052/vacuum.json --hard .52 --execute
~~~

[.52 audited metrics](results/warm-coupled-two-packet-hard052-2026-11-26.json);
[exact .52 terminal state](results/warm-coupled-two-packet-q64-hard052-2026-11-26.json.gz);
[matched terminal comparison](results/warm-coupled-two-packet-matched-hardness-2026-11-26.json).
