# Compare the q16/q32 shifts on matched four- and eight-block clocks

This restarts the actual 12 MeV q16 and q32 prompt-plus-elastic
and signed pair histories to the same 10.442 MeV endpoint.
Earlier unknown exterior **counts** are carried once in the
number audit, not inserted as reconstructed occupations or
gravity. Newly emitted exterior cohorts are known and are
carried through the restart. The method and initial states are
the same as in the [four-block support control](INDEPENDENT-WARM-EARLY-SERIAL-SUPPORT-2026-11-26.md).

At each common physical endpoint let `D_n = log(Ω(q32,n)/Ω(q16,n))`.
The difference `D8-D4` is exactly the difference of the two
within-support time shifts. This tests partial cancellation in
the support comparison; an absolute time-grid shift larger
than a support difference does not by itself measure their
interaction. The same doubled grid has to be run for both
supports before making the comparison.

The separate shape audit computes count, kinetic energy, RMS
momentum, fourth moment, pressure and CDFs of the **retained**
states on independent positive quadratures. It does not
infer momentum moments from missing old tails or use the late
inclusive rest-count proxy as a shape normalization. The q16
32-block state is a separate time-refinement comparator, not
matched against eight-block q32 as a support result.

## Matched result

| At visible T = 10.442 MeV | q16, 4 | q32, 4 | q16, 8 | q32, 8 |
| --- | ---: | ---: | ---: | ---: |
| Inclusive scalar fixed-count Ωh² | .0001622275220 | .0001622691625 | .0001624413069 | .0001624880750 |
| Inclusive stable fixed-count Ωh² | 84.36601067 | 84.35632692 | 84.31629325 | 84.30541703 |

At 8 blocks, the q32 correction is about **+0.02879% scalar**
and **−0.01290% stable**. The signs agree with 4 blocks at
all four probes; the interaction is 12–15% of the 4-block
support log-shift. Thus the larger absolute changes with
time refinement cancel substantially in this *matched*
support comparison. This is not a support- or time-continuum
proof, nor a statement about another visible temperature.

In the retained shapes, q32 relative to q16 at 8 blocks gives
approximately +0.0063% in scalar RMS momentum but +0.56% in
the scalar fourth moment. The analogous stable changes are
+0.117% and +0.88%. The maximum CDF shifts are about
0.00083 (scalar) and 0.00069 (stable). The scalar RMS support
correction is notably smaller at 8 blocks than at 4, so the
count-level partial cancellation must not be generalized to
all shape statistics. The archived moments are separate from
any measured transfer function or future production prediction.

## Reproduction

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
S16="$D"/50MeV-from-empty-pair-quantum-shape-2026-11-26.json
S32="$D"/50MeV-from-empty-pair-quantum-q32-shape-2026-11-26.json
for x in 16 32; do
  if [ "$x" = 16 ]; then seed="$S16"; else seed="$S32"; fi
  python -W error warm_coupled_serial_restart_reference.py EOS.dat "$seed" "q$x-t8" \
    --stop 1.15 --steps 8 --local-prompt-jacobian
done
python -W error warm_early_support_time_interaction_reference.py \
  "$D"/50MeV-from-empty-pair-quantum-q32-audit-2026-11-26.json "$S16" "$S32" \
  "$D"/warm-early-q{16,32}-serial-checkpoint-2026-11-26.json q16-t8/checkpoint.json q32-t8/checkpoint.json
python -W error warm_early_retained_shape_interaction_reference.py \
  "$D"/warm-early-q{16,32}-serial-checkpoint-2026-11-26.json q16-t8/checkpoint.json q32-t8/checkpoint.json \
  "$D"/warm-early-q16-time32-checkpoint-2026-11-26.json.gz
~~~
