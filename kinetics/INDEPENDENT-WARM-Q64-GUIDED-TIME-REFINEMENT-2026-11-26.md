# Refine time on the same evolving q64 support through 10.137 MeV

The [two-block q64 continuation](INDEPENDENT-WARM-Q64-GUIDED-SERIAL-2026-11-26.md)
did not quantify time error on the wider grid. Start all three
comparisons from the same accepted q64 snapshot at 10.238 MeV.
Keep the same physical endpoint, grid, boundary handling,
elastic panels, pair kicks and prompt equations, and split
each archived block into one, two or four substeps.

At intermediate times, the initial Newton coordinates come
from a linear interpolation of the pre-existing narrower
entropy coefficients. They are not stored as measured physical
states and do not enter the accepted q64 right-hand side, moment
targets or ledgers. Archived q32 endpoints can still be used
directly as the initial guesses where they match. Each of the
four- and eight-block solves accepts after three Newton updates.

The 2→4 common-support number (L^1) changes are
(4.88	imes10^{-4}) (scalar) and (2.98	imes10^{-4})
(light), relative to the respective terminal two-block number.
The 4→8 changes are (2.47	imes10^{-4}) and (1.50	imes10^{-4}).
The ratios for count and energy all lie between 1.97 and 1.99,
consistent with first-order time behavior on this finite
grid and interval. They are much larger than the separately
observed q32-to-q64 support changes on the original two blocks.
That does not mean support error is bounded on every time grid
or at much later temperatures.

Only the selected coupled operators and the same known
collisionless gravity are included; old exterior distributions
and additional physics are still not reconstructed.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Q64=results/warm-extended-coupled-q64-3-blocks-2026-11-26.json.gz
Q32=results/warm-light-occupation-next-q32-tape-2026-11-26.json.gz
OBS=results/warm-observed-suboperator-q32-tape-2026-11-26.json.gz
BASE=results/warm-guided-q64-next-tape-2026-11-26.json.gz
python -W error warm_q64_guided_refinement_reference.py EOS.dat "$Q64" "$Q32" "$OBS" "$BASE" run4 --splits 2
python -W error warm_q64_guided_refinement_reference.py EOS.dat "$Q64" "$Q32" "$OBS" "$BASE" run8 --splits 4
python -W error warm_q64_ladder_endpoint_audit_reference.py "$BASE" \
 run4/run/step-4.json run8/run/step-8.json
~~~

[Independent endpoint audit](results/warm-guided-q64-ladder-audit-2026-11-26.json);
[four-block endpoint](results/warm-guided-q64-4-blocks-endpoint-2026-11-26.json.gz)
and [summary](results/warm-guided-q64-4-blocks-summary-2026-11-26.json);
[eight-block endpoint](results/warm-guided-q64-8-blocks-endpoint-2026-11-26.json.gz)
and [summary](results/warm-guided-q64-8-blocks-summary-2026-11-26.json).
