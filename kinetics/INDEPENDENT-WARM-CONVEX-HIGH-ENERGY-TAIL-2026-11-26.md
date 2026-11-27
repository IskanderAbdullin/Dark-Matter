# Track creation of energy above a high threshold, not just passage through a heavy parent

For a fixed physical energy cut C, define w(E) = max(E−C,0).
The energy moment W is the sum of w(E) over **all** species.
For a binary merger the per-event increase is exactly

~~~text
d(x,y) = min(x,y,C,max(x+y−C,0)).
~~~

The reverse decay subtracts the same amount. A two-to-two
collision changes W by the difference of this defect between
the two sides, at fixed total energy. Without collisions,
redshift can only decrease the comoving W for a fixed physical
cut. This is a distinct observable from the gross energy entering
a heavy parent each time it forms. For the earlier [synthetic
hard packet control](INDEPENDENT-WARM-HIGH-INVERSE-FROZEN-PRESSURE-TEST-2026-11-26.md),
the 100 GeV merger adds only 4.8e-7 of its gross parent energy
to W at C=0.71482 GeV: its incoming light particle already
carried almost all the high energy. This observation is not
a bound on how often the event recurs.

## Resolved frozen rates at 7.394049 MeV

Evaluate this moment on the final corrected q64-48 state.
For prompt decay/inverse, use exact angles in both directions
and a stationary heavy occupation. Compute the signed rate two
ways: from the merger defect and from the emitted/consumed
daughter W, since steady-parent production cancels at each parent
energy. For the common scalar-pair chart and retained below-pole
elastic operator, accumulate the event weights before projecting
them into moments. Split angular/energy integration at the
observable cuts and compare three orders. The optional diagnostic
is not used for the accepted time integration.

| C (GeV) | Prompt with vacuum incoming exterior | Common scalar pair | Retained below-pole elastic |
| ---: | ---: | ---: | ---: |
| .10 | 4.888 | 1.193e-5 | 6.107e-6 |
| .30 | 1.255e-11 | 1.483e-17 | 9.899e-19 |
| .50 | 4.658e-22 | 7.748e-30 | 0 |
| .71482 | 2.008e-36 | 4.035e-44 | 0 |

The table entries are signed W production rates divided by
retained dark energy times H. The largest cut corresponds to
a heavy momentum cut .70779 GeV. At the last two cuts all
represented daughter energies are below the cut, so closed
stationary prompt has zero signed W exactly. The open prompt
result there is the independently counted positive excess
of newly outgoing daughters, not a subtraction of near-equal
large numbers. At .71482 GeV the q64 gross inverse parent
energy rate remains 4.117e-32 in the same units. No high
outgoing elastic states are allowed by the closed elastic model;
the zero is geometric, not a physical high elastic rate limit.
The uppermost pair result moves about 5% between the two finer
orders, so read it at the order-of-magnitude level.

## A bound on the already stored outside energy at the same cut

For a positive collisionless birth cohort with scale ratio r,
mass m and unknown birth energy e, the present energy is
F(e)=sqrt(m²(1−r²)+r²e²). Both F and max(F−C,0) are convex.
Use the known mean and birth support to bound the latter:
Jensen below, the chord above on a finite interval, or the
upper slope r above the floor on an unbounded high interval.

For the actual energy above C, majorize F(e) 1{F(e)≥C} by
a concave piecewise-linear envelope. Below the threshold
preimage, join the supporting floor to the point with energy C;
above it use slope r. If the entire support already starts
above the cut, use slope r from that floor; if it lies below the
cut, charge zero. This uses only cohort moments, not a chosen
intra-cohort shape. Random finite and semi-infinite controls
check both envelopes.

At .71482 GeV the q64-48 stored receipts can hold at most
9.108e-13 of the conditional late full energy allowance above
the cut, counting both light and scalar. The q32-48 receipts
give 3.028e-12, and the Jensen lower in both cases is zero.
These are bounds for **the recorded collisionless cohorts only**.
They do not limit an unarchived incoming population, future
reacting re-entry or repeated growth. The frozen-rate table
and the receipt bounds cannot be read as a coupled tail bound.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/warm-high-inverse-frozen-pressure-test-2026-11-26.json
a=results/warm-corrected-q32-time48-step48-2026-11-26.json.gz
b=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
python -W error warm_prompt_excess_energy_tail_reference.py EOS.dat "$p" "$a" "$b"
python -W error warm_pair_chart_excess_energy_tail_reference.py EOS.dat "$p" "$b"
python -W error warm_elastic_excess_energy_tail_reference.py EOS.dat "$p" "$b"
python -W error warm_exterior_convex_tail_budget_reference.py \
  results/warm-corrected-anchor-signed-floor-2026-11-26.json "$a" "$b"
~~~

[Prompt](results/warm-convex-high-tail-prompt-2026-11-26.json);
[pair](results/warm-convex-high-tail-pair-2026-11-26.json);
[elastic](results/warm-convex-high-tail-elastic-2026-11-26.json);
[stored cohort bounds](results/warm-convex-high-tail-exterior-2026-11-26.json).
