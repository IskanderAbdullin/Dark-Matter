# Separate small pair feedback from the frozen scheduling error

The [frozen operator schedules](INDEPENDENT-WARM-FROZEN-SPLIT-ORDER-2026-11-26.md)
produce different retained-plus-outside counts at the same
normalization and physical shell geometry. Here compute the
leading count response to the explicit pair operator without
using either split schedule for that derivative.

Introduce a **numerical** multiplier λ of the whole signed pair
operator. Evolve only the prompt-plus-elastic collision map at
λ=0 with the L-stable integrator on the identical frozen warm
geometry. At every archived point evaluate the signed pair
rate on that collision-only state. Since the fast operator
conserves each species number, the first derivative of the
number change with respect to λ at zero is simply the time
integral of those rates. No derivative of unknown exterior
occupations is needed for that leading calculation.

Two positive collision-only resolutions (22 and 43 points,
clustered near the start) and the two integration rules give
relative spread 5.2×10⁻⁶ for the two number changes. The signed
first variation, divided by the respective initial counts,
is about **+0.00347356 scalar**, **−0.00139073 stable** for the
fixed normalized interval of 0.006. The scalar net rate per
initial scalar count changes from 0.5511 to 0.5802 per unit
log-expansion interval. The inverse/forward ratio changes
from 0.8399 to 0.8315. The full prompt-plus-elastic reaction
alone changes neither count.

Divide the measured finite-strength count difference by λ
to compare with that derivative. At λ=.03:

| Schedule | Relative shortfall vs leading count derivative, 1 substep | At 4 substeps |
| --- | ---: | ---: |
| Pair half, collision, pair half | 2.368% | 0.636% |
| Collision half, pair, collision half | 0.045% | 0.049% |

The λ=.10 case and the archived λ=1 case are also recorded.
Two outcomes matter. A scheduling discrepancy is still present
when feedback from pairs is small. And the near agreement of
the reverse order at small λ is a control of **this frozen
leading count response**, not evidence that it is always a
more accurate or cheaper whole-trajectory solver. The newer
collision stepper alone does not settle the issue. Both tests
leave the actual default warm calculation unchanged. This
does not vary a measured coupling or estimate a physical
uncertainty in an expanding history.

~~~sh
export OPENBLAS_NUM_THREADS=1
S=results/50MeV-from-empty-pair-quantum-shape-2026-11-26.json
python -W error warm_frozen_leading_pair_reference.py EOS.dat "$S" >leading.json
for x in .1 .03; do
  python -W error warm_frozen_split_order_reference.py EOS.dat "$S" --method sdirk2 \
    --schemes pair-collision-pair collision-pair-collision --substeps 1 4 --pair-strength "$x" >"split$x.json"
done
python -W error warm_frozen_leading_pair_audit_reference.py EOS.dat "$S" leading.json \
  results/warm-frozen-split-sdirk-2026-11-26.json split.1.json split.03.json
~~~
