# Complete the matched warm 4, 8, 16 and 32-block comparison

The 50 MeV warm q16 and q32 continuations now both reach
relative expansion b=1.15 and 10.442 MeV at each resolution.
The latest q32 run reused the verified atomic checkpoint in
blocks of eight without changing its plan. Count earlier
exterior inventory once analytically; the older momentum
distributions were not restored to the evolution. Newly
tracked exterior cohorts are carried across restarts.

| Matched 32-block endpoint | Inclusive scalar Ωh² proxy | Inclusive stable Ωh² proxy |
| --- | ---: | ---: |
| q16 | .00016258786 | 84.28221086 |
| q32 | .00016264020 | 84.27003819 |

Let D be ln(q32/q16) at a fixed probe and time resolution.
On all 16 matched 16→32 probes the sign stays positive for
scalar and negative for stable. At the end, D32 is
(+3.21884×10⁻⁴, −1.44438×10⁻⁴). Its increment from D16 is
3.77% and 3.81% of the respective finest contrasts.
The increment shrinks by a factor near .553 from the
8→16 step; the apparent endpoint order changes from about
0.51 to 0.85. Do not infer a continuum order from this.

Retained shapes have separate sensitivity. At the same
32-block endpoint the scalar RMS-momentum support contrast
crosses zero, from +5.64×10⁻⁶ at 16 to −2.44×10⁻⁵ at 32,
while the scalar fourth moment still shifts by .542%.
The stable fourth moment shifts by .890%. The maximum
absolute retained CDF support difference at the ten checked
thresholds is .000968 for scalar and .000745 for stable.
These are early retained-only moment checks, not evidence
for full momentum convergence or a late relic. The warm
stable proxy remains far from .120.

The checked-in gzip file is the entire atomic endpoint.
Run the count and shape comparators with `--base-blocks 16`
and the ladder with `--finest`; earlier outputs are unchanged.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
A="$D"/warm-early-q16-time16-checkpoint-2026-11-26.json.gz
B="$D"/warm-early-q32-time16-checkpoint-2026-11-26.json.gz
C="$D"/warm-early-q16-time32-checkpoint-2026-11-26.json.gz
E="$D"/warm-early-q32-time32-checkpoint-2026-11-26.json.gz
python -W error warm_early_support_time_interaction_reference.py \
 "$D"/50MeV-from-empty-pair-quantum-q32-audit-2026-11-26.json \
 "$D"/50MeV-from-empty-pair-quantum{,-q32}-shape-2026-11-26.json \
 "$A" "$B" "$C" "$E" --base-blocks 16 > time1632.json
python -W error warm_early_retained_shape_interaction_reference.py \
 "$A" "$B" "$C" "$E" "$C" --base-blocks 16
python -W error warm_early_support_ladder_summary_reference.py \
 "$D"/warm-early-support-time-interaction-2026-11-26.json \
 "$D"/warm-early-support-time-8-16-interaction-2026-11-26.json --finest time1632.json
~~~
