# Continue the guided q64 solution through 10.137 MeV

Apply the [verified borrowed-initialization idea](INDEPENDENT-WARM-Q64-TRANSPLANTED-NEWTON-2026-11-26.md)
to two new q64 blocks beyond the last accepted q64 state at
10.238 MeV. The already archived q32 endpoints supply only
the starting coordinates for two nonlinear solves on wider
q64 support. The accepted targets, pair kicks, elastic panels,
prompt dynamics, redshift and old collisionless gravity stay
on q64. If the borrowed solve fails, the wrapper can fall back
to the existing wider solve without accepting the trial guess.

Both new nonlinear solves accept after four Newton updates,
without fallback. Their accepted normalized residuals are
(4.62	imes10^{-12}) and (3.03	imes10^{-12}). The
independently checked retained-plus-new-receipt number defect
is below (10^{-14}), and inherited cohorts rebase exactly
on both checkpoints. The output branches are recorded with
their source and guide digests, not just final summaries.

The terminal visible temperature is 10.136743 MeV. On the
common shells the q32/q64 count (L^1) differences are
(7.21	imes10^{-7}) for scalars and (1.29	imes10^{-7})
for light fermions, relative to the corresponding q32 count.
The explicitly admitted q64-only region has resident scalar
and light fractions (1.10	imes10^{-11}) and
(2.18	imes10^{-11}). The largest retained light occupation
in these two **reconstructed endpoint** snapshots remains
below (5.79	imes10^{-5}). No pointwise statement about the
physical distribution within unknown cells or between endpoints
is implied.

This is still the selected coupled operator set, with no
new physical source, actual heavy phase, or older untracked
exterior reconstructed. The number of nonlinear iterations
is a measured solver result, not a timing or general convergence
guarantee.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Q64=results/warm-extended-coupled-q64-3-blocks-2026-11-26.json.gz
GUIDE=results/warm-observed-suboperator-q32-tape-2026-11-26.json.gz
python -W error warm_q64_guided_serial_reference.py EOS.dat "$Q64" "$GUIDE" guide-run
python -W error warm_q64_guided_restart_audit_reference.py "$Q64" "$GUIDE" \
 results/warm-guided-q64-next-tape-2026-11-26.json.gz
~~~

[Exact atomic tape](results/warm-guided-q64-next-tape-2026-11-26.json.gz),
[solver summary](results/warm-guided-q64-next-summary-2026-11-26.json)
and [independent restart audit](results/warm-guided-q64-next-audit-2026-11-26.json).
