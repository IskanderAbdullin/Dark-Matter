# An equal initial tail can still grow or shrink

The [same-center width comparison](INDEPENDENT-WARM-LATE-NEAR-CUT-SHAPES-2026-11-26.md)
holds the archive moments but changes initial W. To test W
as an additional predictor, solve for the hard-packet center
in two other width families such that the total W is the same.
Both keep the original exterior birth count and energy of
each recorded cohort and the 1.05 GeV closed pad at the same
visible temperature. Measured initial species number and
energy coincide exactly at stored precision, and total
W/rho matches 2.0e−16 to 9e−14 relative.

| Hard center (GeV) | Half-width (GeV) | Signed current divided by rho H |
| --- | --- | ---: |
| .7134635061 | .00250 | +7.05153e−10 |
| .7148726662 | .00025 | −6.98683e−10 |

All three independent angular quadrature orders agree
on the signs. Conserving the same particle counts, the
same total energy, and the same total convex excess
does not fix that excess's initial rate in this closure.
The individual species' W and detailed momentum populations
differ; this does not test a closure retaining those individual
W moments, an arbitrary physical exterior, or a full reaction network.

## Lifetime versus time grid on the separate .713 GeV pair

Repeat the earlier .713 GeV broad and narrow reconstructions
with a time-dependent parent initialized at its stationary
occupation. Retain the 7,730 parent nodes and use one backward
Euler step at each of the six previously saved endpoints.
To isolate the parent lifetime effect, compare against
stationary runs with **one** step per interval, not only
against the earlier three-substep runs. Per-species total
daughter-plus-parent numbers and total massive energy pass
the conservation checks, and the heavy equation passes a
roundoff-scaled independent residual check.

The effect of finite parent lifetime on W at the matched
clock is at most 3.59e−8 and 3.19e−7 of initial W for the
5 MeV and 2.5 MeV half-width cases respectively. The signs
of the sampled currents agree on both clocks and closures:
the broad case loses W throughout, the narrow one initially
gains it and the saved terminal current is negative. The
one-versus-three stationary clock change in terminal W
is +2.78% (broad) and −1.24% (narrow) relative to the
three-step trajectory. Those percent-level differences
should not be assigned to the parent lifetime. The parent
starts in equilibrium; a different parent population or an
expanding, continuing source is not tested by this comparison.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
a=results/warm-corrected-anchor-signed-floor-2026-11-26.json
b=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
python -W error warm_late_equal_excess_shape_rate_reference.py EOS.dat "$a" "$b"
# Extract 713-3 and 713-narrow-3 as in the near-cut count note.
for case in 713 713-narrow; do
  if [ "$case" = 713 ]; then w=.005; else w=.0025; fi
  python -W error warm_late_near_crossing_shape_evolution_reference.py EOS.dat "$a" "$b" "$case-1" --hard .713 --halfwidth "$w" --substeps 1
  python -W error warm_late_finite_parent_prompt_reference.py EOS.dat "$case-parent" "$case-3"/step-{0,1,2,3,4,5,6}.json --initial equilibrium
done
python -W error warm_late_near_cut_finite_parent_audit_reference.py --matched 713-1 713-narrow-1 --finite 713-parent 713-narrow-parent --fine 713-3 713-narrow-3
~~~

[Equal-W phases and frozen currents](results/warm-late-equal-excess-shape-rate-2026-11-26.json);
[matched parent and clock audit](results/warm-late-near-cut-finite-parent-match-2026-11-26.json).
