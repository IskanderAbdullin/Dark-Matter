# A q64 history-based guess does not need the next q32 shape

The [q32-guided solver](INDEPENDENT-WARM-Q64-GUIDED-SERIAL-2026-11-26.md)
required an existing narrower shape to initialize each wider
block. Test an alternative for the next interval, from
10.136743 to 10.087326 MeV. Extrapolate entropy coefficients
from the last two already accepted **q64** states at the
scheduled log-expansion interval. The extrapolation is
only an initial guess for the full nonlinear implicit q64
solve. Do not extrapolate or replace any of the accepted
physical moment targets or force columns.

Compare to an independent replay initialized by the archived
q32 endpoint. That snapshot supplies the common time
interval to both tests but **only** the second solve
receives its distribution as a guess. The q64-history
solve accepts in three Newton updates to (5.33	imes10^{-12})
normalized residual; the q32-guided replay accepts in four.
The retained q64 number and energy moments of the two
solutions agree to under (3.7	imes10^{-15}) in the
tested relative (L^1) metrics.

This is one measured continuation block, with the original
prompt, selected elastic, pair and exterior-gravity operators.
It does not prove a universal stable predictor, an elapsed-time
speedup, or physical validity of the selected model at later
times. Keep accepted physical checkpoints for continuation.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_q64_history_predictor_reference.py \
 EOS.dat results/warm-guided-q64-next-tape-2026-11-26.json.gz \
 results/warm-q32-lower-guide-first-block-2026-11-26.json.gz run-history
~~~

[Diagnostics](results/warm-q64-history-predictor-2026-11-26.json),
[history-based endpoint](results/warm-q64-history-predictor-endpoint-2026-11-26.json.gz)
and [independently guided endpoint](results/warm-q64-history-predictor-guided-control-2026-11-26.json.gz).
