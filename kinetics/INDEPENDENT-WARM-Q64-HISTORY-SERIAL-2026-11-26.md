# Continue independently on q64 to 9.748 MeV using only its own history

The [single-block history predictor control](INDEPENDENT-WARM-Q64-HISTORY-PREDICTOR-2026-11-26.md)
showed that the wider solver need not use a future narrow
distribution. Start from the accepted q64 state at 10.087326 MeV
and its predecessor, fix the next seven log-expansion intervals,
and use each pair of accepted q64 shapes only to seed the next
nonlinear solve. New accepted targets are computed by the
original wider physical operators. If a history trial fails,
only the accepted wider fallback may supply a checkpoint.

Every block accepted with the history-based initialization in
three Newton updates, with no fallback. The endpoint visible
temperature is 9.747990 MeV; the relative expansion for these
seven blocks is (1.04^{7/8}). An independent restart audit
checks the parent hash of every accepted snapshot and the
stored source digest, rebases the original exterior cohorts,
and verifies the cumulative retained-plus-new-receipt number
to (2	imes10^{-14}). The original local energy and work
ledgers are recorded alongside the state. No narrower future
distribution is required in this serial driver.

The accepted wider snapshots include both selected pair
kicks, prompt transfer, selected elastic scattering and
the existing gravity treatment. They do not include new
physical visible injection, an evolved heavy phase, or late
cosmology. A measured iteration count is not a runtime or
general convergence bound.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
zcat results/warm-q64-history-serial-seed-2026-11-26.json.gz > seed.json
STEP=$(python -c 'print(1.04**(1./8))')
python -W error warm_q64_history_serial_reference.py EOS.dat seed.json history-run \
 --step-expansion "$STEP" --blocks 7
python -W error warm_q64_history_serial_audit_reference.py EOS.dat seed.json \
 results/warm-q64-history-serial-a-2026-11-26.json.gz \
 results/warm-q64-history-serial-b-2026-11-26.json.gz
~~~

[Seed](results/warm-q64-history-serial-seed-2026-11-26.json.gz),
[blocks 1–3](results/warm-q64-history-serial-a-2026-11-26.json.gz),
[blocks 4–7](results/warm-q64-history-serial-b-2026-11-26.json.gz),
[solver summary](results/warm-q64-history-serial-summary-2026-11-26.json)
and [independent audit](results/warm-q64-history-serial-audit-2026-11-26.json).

## Independent matched q32 calculation

A separate original-solver q32 calculation also reaches
9.747990 MeV on the same absolute grid in log expansion.
Each corresponding checkpoint temperature and shared
momentum boundary matches the independent q64 trajectory.
At the terminal point the common-shell count (L^1)
differences are (1.20	imes10^{-6}) (scalar) and
(2.44	imes10^{-7}) (light), relative to the respective
retained q32 counts. The extra q64 region holds
(1.03	imes10^{-11}) and (2.13	imes10^{-11}) of its
represented counts. The q32 calculation never supplies
history-seed shapes to these q64 blocks.

Add `--narrow` followed by both linked q32 tapes to reproduce
the [matched audit](results/warm-q64-history-serial-matched-2026-11-26.json).
The q32 [first four](results/warm-q32-lower-matched-a-2026-11-26.json.gz)
and [last four](results/warm-q32-lower-matched-b-2026-11-26.json.gz)
snapshots and [exact terminal checkpoint](results/warm-q32-lower-endpoint-2026-11-26.json.gz)
are retained. This checks support on the same selected time
grid, not continuum convergence at the lower endpoint.
