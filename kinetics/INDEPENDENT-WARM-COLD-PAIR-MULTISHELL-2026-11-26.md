# Lift the reverse pair events into six or nine cold shells

A [single-cell partition](INDEPENDENT-WARM-COLD-PAIR-PARTITION-2026-11-26.md)
removes the cold-cell spike but leaves the same incoming-scalar projection
in nearby cells. Extend the disjoint event partition to the union of the
first six or nine light cells. Split the tagged radial integral at every
shell edge and project each event to all four legs. If both light products
are below the common cut, divide by two, even when they land in different
shells. The quadrature on the complement and the entire forward rate stay
unchanged. The recovered total gain per lifted cell is cross-checked
against tagging only that cell; this independently tests copy symmetry,
not just conservation of aggregate number or energy.

The five frozen tests use actual q32 and q64 checkpoints spanning
8.133–7.724 MeV. With the nine-cell cut, the sampled momentum range ends
near 1.27–1.33 MeV, and the lift processes about 0.025–0.027% of all
reverse events. Near the spike the first four native pair contributions
are positive on both grids, about +9.0e-5, +9.4e-5, +5.3e-5 and +2.6e-5
per Hubble time. The nine-cell rule instead gives negative rates near
−1.1e-5 in all four. The two tagged quadratures agree within 1.0e-7
fractionally across the corrected signed cell rates on these snapshots.
Their variation in the complete signed-rate shell L1, normalized to
the native per-species per-moment L1, is below 8.2e-8. The total collision
rate still depends on the unrefined complement. An experimental
twelve-cell cut did not pass the copy-symmetry tolerance and is not
a supported setting here.

Repeat the independent q32 replay from identical actual predecessors.
At 7.754805 MeV the original first-cell mean slope is +1.43547e-4
per ln(a); the one-, six- and nine-cell variants each give
approximately −4.7185e-5. The largest normalized whole-grid change
in number, energy and pressure over this same interval is 3.4e-10,
2.2e-8 and 6.7e-8 across the respective partition sizes, taking the
largest of six moment/species comparisons for each variant.
The nearest neighboring cells also lose the artificial gains.
The truncated kinetic network and the evolution time grid are unchanged.
No final relic prediction follows from this comparison.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# Prepare the native replay inputs as in the single-cell test.
for n in 6 9; do
  python -W error warm_cold_pair_spike_replay_reference.py EOS.dat replay replay/32-cut-$n --q 32 --cells "$n"
done
python -W error warm_cold_pair_multishell_audit_reference.py EOS.dat \
  results/warm-recoil-crossing-q{32,64}-terminal-2026-11-26.json.gz \
  results/warm-cold-pair-spike-q{32,64}-step12-2026-11-26.json.gz \
  results/warm-cold-pair-spike-q64-step13-2026-11-26.json.gz
~~~

[Frozen rate comparisons](results/warm-cold-pair-multishell-audit-2026-11-26.json.gz);
[actual replay results](results/warm-cold-pair-multishell-spike-q32-2026-11-26.json);
[exact nine-cell replay output](results/warm-cold-pair-multishell-q32-9-lifted-2026-11-26.json.gz).
