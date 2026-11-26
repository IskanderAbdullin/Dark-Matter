# Carry the matched q16/q32 comparison through sixteen blocks

Extend the [complete eight-block comparison](INDEPENDENT-WARM-EARLY-SUPPORT-TIME-INTERACTION-2026-11-26.md)
using the same q16 and q32 12 MeV restart spectra, local
collision algorithm, and visible radiation clock. Both runs
complete the prescribed expansion factor 1.15 at a visible
temperature of 10.442 MeV. Earlier unresolved source tails
are included once as counts only; their unknown occupations
are not reconstructed. Newly projected positive exterior
cohorts are tracked in the restarted trajectories.

| At the endpoint | Eight blocks, q16 | Eight blocks, q32 | Sixteen, q16 | Sixteen, q32 |
| --- | ---: | ---: | ---: | ---: |
| Inclusive scalar fixed-count Ωh² | .0001624413069 | .0001624880750 | .0001625401869 | .0001625905439 |
| Inclusive stable fixed-count Ωh² | 84.31629325 | 84.30541703 | 84.29329792 | 84.28158707 |

At 16 blocks the widened support changes the inclusive count
by **+0.0310% scalar** and **−0.0139% stable** at this point.
The sign agrees with the 8-block comparison at all eight
matched temperatures. The change in the support log contrast
from 8 to 16 blocks is smaller than each within-support time
shift, and is about 7.6–7.7% of the 8-block support contrast
at the endpoint. Comparing directly with the common four
original endpoints, all three supported contrast signs persist.
At the endpoint the 8→16 contrast difference is about 0.70
of the 4→8 difference (apparent order about 0.51).
This is partial cancellation on finite grids, not a bound
on the continuum answer.

For the retained phase-space profiles the 16-block scalar RMS
support contrast is now approximately (+5.64\times10^{-6})
relative, while the scalar fourth-moment contrast is still
about **+0.55%**. Stable RMS/fourth-moment support changes are
about +0.120%/+0.89%. These are independently evaluated from
positive reconstructed shapes, not synthesized from the final
count. The old outside shape remains unknown. The q16
32-block run is used as its own time-control only.

This compares early warm numerical trajectories, not a final
relic density, a continuation to 3.43 MeV, or an observational
transfer calculation.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
S16="$D"/50MeV-from-empty-pair-quantum-shape-2026-11-26.json
S32="$D"/50MeV-from-empty-pair-quantum-q32-shape-2026-11-26.json
for x in 16 32; do
  if [ "$x" = 16 ]; then seed="$S16"; else seed="$S32"; fi
  python -W error warm_coupled_serial_restart_reference.py EOS.dat "$seed" "q$x-t16" \
    --stop 1.15 --steps 16 --local-prompt-jacobian
done
A="$D"/warm-early-q16-time8-checkpoint-2026-11-26.json.gz
B="$D"/warm-early-q32-time8-checkpoint-2026-11-26.json
python -W error warm_early_support_time_interaction_reference.py \
 "$D"/50MeV-from-empty-pair-quantum-q32-audit-2026-11-26.json "$S16" "$S32" \
 "$A" "$B" q16-t16/checkpoint.json q32-t16/checkpoint.json --base-blocks 8 >time816.json
python -W error warm_early_retained_shape_interaction_reference.py \
 "$A" "$B" q16-t16/checkpoint.json q32-t16/checkpoint.json \
 "$D"/warm-early-q16-time32-checkpoint-2026-11-26.json.gz --base-blocks 8
python -W error warm_early_support_ladder_summary_reference.py \
 "$D"/warm-early-support-time-interaction-2026-11-26.json time816.json
~~~
