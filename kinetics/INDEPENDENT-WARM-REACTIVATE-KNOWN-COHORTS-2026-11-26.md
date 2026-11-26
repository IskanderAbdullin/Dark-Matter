# Allow the known received exterior cohorts to react in an explicit control

The [matched q64 comparison](INDEPENDENT-WARM-MATCHED-EXPLICIT-SUPPORT-COUPLED-2026-11-26.md)
kept the earlier received exterior cohorts as gravity only,
with empty explicit padding. Test one admissible alternative
for those known cohorts, without assigning a distribution to
the older untracked population.

Place each known cohort at the momentum corresponding to
its mean energy at birth and redshift that momentum to the
actual 10.391 MeV anchor. Transfer its number and kinetic
energy into the appropriate q64 shell moments and remove the
same entry from the collisionless receipt ledger. Reconstruct
the transferred shells with the same positive two-moment
closures. This is a one-representative-momentum prescription:
it preserves the stored birth number and energy and corresponds
to the archived Jensen-lower present energy, not the true
unknown spread. Both initial species numbers and modeled
total energy are unchanged by the transfer.

All 132 known cohorts fit in the added support. The numbers
transferred are (8.82	imes10^{-12}) (scalar) and
(1.26	imes10^{-12}) (light) of initial represented species
numbers. Evolve both selected-channel branches over the same
two blocks to (b=1.001). The resulting shared-core number
(L^1) differences are (4.68	imes10^{-11}) and
(4.90	imes10^{-11}), normalized to the initial species
counts. The added-shell resident populations both relax back
close to the comparison with empty initial extra shells.
The native conservation ledgers close below (4	imes10^{-15}).

This tests sensitivity to one phase-space representation of
the recorded cohorts and lets those particles react explicitly.
It is not a bound on arbitrary within-cohort momenta, a result
for inherited untracked tails, or long-term abundance proof.

~~~sh
zcat results/warm-matched-explicit-support-q64-2-blocks-2026-11-26.json.gz > baseline.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_reactivate_known_cohorts_reference.py \
 EOS.dat run-reactivated baseline.json --execute
~~~

[Diagnostics](results/warm-reactivate-known-cohorts-2026-11-26.json);
[exact endpoint](results/warm-reactivate-known-cohorts-q64-checkpoint-2026-11-26.json.gz).
