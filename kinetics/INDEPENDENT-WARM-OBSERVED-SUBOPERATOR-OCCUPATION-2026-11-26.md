# The sampled light peak can rise inside a block while its endpoints fall

Continue both previous [matched saved grids](INDEPENDENT-WARM-LIGHT-OCCUPATION-NEXT-MATCHED-STAGE-2026-11-26.md)
from their 10.238 MeV endpoints for two further matched
blocks through another 1% of expansion, to 10.137 MeV.
Attach read-only observers at redshift, the pair
reconstructions, and the accepted output of the implicit
collision update. Store entropy-array digests and
occupation extrema. The explicit pair midpoint predictors
are labeled as trials, not accepted physical states.
A separately executed uninstrumented q16 first block
has an atomic checkpoint identical byte for byte.

| q | Block | Entrance light peak | Largest observed nontrial peak | Exit light peak |
| --- | ---: | ---: | ---: | ---: |
| 16 | 1 | 5.954630×10⁻⁵ | 5.956387×10⁻⁵ | 5.936003×10⁻⁵ |
| 16 | 2 | 5.936003×10⁻⁵ | 5.937363×10⁻⁵ | 5.915779×10⁻⁵ |
| 32 | 1 | 5.795670×10⁻⁵ | 5.803189×10⁻⁵ | 5.786169×10⁻⁵ |
| 32 | 2 | 5.786169×10⁻⁵ | 5.794213×10⁻⁵ | 5.784931×10⁻⁵ |

The largest value in each row follows the **accepted
first pair half-kick**. The implicit collision then reduces
the peak. Thus an endpoint-only decreasing sequence is
not a maximum envelope even in the same discrete method.
The largest q32 within-block excess is (8.044×10^{-8})
in absolute occupation. It is tiny compared with the
analytic half-per-spin criterion, but these checks do not
bound unsampled internal nonlinear states, interpolations,
earlier omitted exterior occupations or the full history.
Nor does one observed split trajectory distinguish a
physical transient from the ordering of operators.

To reproduce the tapes, wrap the last previous snapshot
in a one-entry snapshots source using the serial driver's
canonical bytes_of. Run the observer on the same directory
twice with defaults (stop 1.01 and two blocks), saving
both returned atomic checkpoints and observation records.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_observed_suboperator_audit_reference.py \
 "$D"/warm-observed-suboperator-q{16,32}-tape-2026-11-26.json.gz \
 "$D"/warm-light-occupation-next-q{16,32}-tape-2026-11-26.json.gz
~~~
