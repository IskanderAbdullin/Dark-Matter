# The activated .58 GeV hard tail relaxes in the selected coupled warm run

Start from the same actual 10.391 MeV q32 shape as the short
[matched support control](INDEPENDENT-WARM-MATCHED-EXPLICIT-SUPPORT-COUPLED-2026-11-26.md).
On the expanded q64 domain, explicitly represent each of 132 previously
collisionless stored exterior cohorts by two moment-matched packets.
Put one momentum at .58 GeV at this anchor and select the other below the
stored mean while respecting its support floor. Match each cohort's count
and mean **birth** energy; reconstruct the resulting shell counts and
lab-frame kinetic moments. The activated receipts no longer appear in the
old collisionless gravity ledger. New outer pair births are still recorded.

Evolve the selected prompt, below-pole elastic and signed scalar-pair
network in two identical expansion intervals matching the vacuum-padded
q64 comparison. Both curves reach the same 10.380162 MeV endpoint at
relative expansion 1.001. Atomic restart hashes, the absolute matched
mesh and clock, and the count including new births are checked.
The final inclusive number defect is 6.4e-15; the two individual energy
ledgers are below 2e-15. The first nonlinear interval accepted five
microsteps and rejected nine; the second accepted one.

The normalized bare inverse energy rate for the heavy comoving band
starting at 1 GeV falls from 6.6391e-21 to 3.7792e-32
(in units of retained dark energy times Hubble), a drop of
eleven orders of magnitude for this realization. The end rate is
within 1.7e-5 relative of the matched vacuum-padded run, but the
shared-shell scalar/light number L1 differences are 1.32e-5 and
4.44e-6. Tail relaxation is not proof the complete histories coincide.
There is no coupled time refinement yet for this realization.

This is a deliberately hard reconstruction of stored receipts, not a
shape-independent bound. It extends the
[isolated prompt-only stress](INDEPENDENT-WARM-DIRECT-PROMPT-STRESS-2026-11-26.md)
to a genuine short run of the selected coupled network. It still does
not evolve an independent heavy spectrum, old unknown exterior phase,
or unselected interactions, and uses the archived native pair
discretization. The output below applies only to the .58 GeV case.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p stress/q64
zcat results/warm-matched-explicit-support-q64-2-blocks-2026-11-26.json.gz > stress/vacuum.json
zcat results/warm-coupled-two-packet-q64-hard058-2026-11-26.json.gz > stress/q64/checkpoint.json
python -W error warm_coupled_two_packet_completed_reference.py EOS.dat stress stress/vacuum.json
# Use a separate empty work directory to evolve from scratch:
python -W error warm_cohort_two_packet_stress_reference.py EOS.dat rerun stress/vacuum.json --hard .58 --execute
~~~

[Audited metrics](results/warm-coupled-two-packet-hard058-2026-11-26.json);
[exact accepted terminal checkpoint](results/warm-coupled-two-packet-q64-hard058-2026-11-26.json.gz).
