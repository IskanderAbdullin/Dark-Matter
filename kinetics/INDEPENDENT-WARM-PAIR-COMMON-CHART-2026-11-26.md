# Put pair gain and loss on the same fermion quadrature

The [cold reverse lift](INDEPENDENT-WARM-COLD-PAIR-MULTISHELL-2026-11-26.md)
uses the fermion chart only for a low-momentum subset. A more direct
alternative is to use the common invariant event measure for both
directions over the whole scalar-pair process. Integrate over an
ordered opposite-charge fermion pair and the two-body scalar final
angles. On that one quadrature, the direct weights are

- inverse gain: (1−f₁)(1−f₂) B₁ B₂;
- forward loss: f₁ f₂ (1+B₁)(1+B₂).

The two spin states on each fermion and the identical-scalar factor
appear for both directions on this chart. There is no subtraction
of tiny boson occupations from unity. A kinematic total-energy bound,
E_L,max = 2 E_phi,max − m_L, limits additional virtual fermion
points for inverse events. These points have zero occupation
and do not add any source phase population. Contributions outside
the material grid go into the original exterior moment ledger.
Every event projects to all four legs, with the proper signs.

The synthetic common-temperature check compares more than 318,000
nonzero interior integrands. The relative difference is below
3.56e-14, and the largest of the first seven net/loss ratios is
5.50e-16 already at the coarser paired quadrature. The full
projected count/energy ledger closes to 8e-17. No exact thermal
stationarity of the artificially truncated exterior is assumed.

On the actual q32 replay to 7.754805 MeV, the first-shell mean
occupation changes at +1.43547e-4 per ln(a) for the original chart,
−4.71851e-5 for the nine-cell lift and −4.71850e-5 for the complete
common chart. The whole-state species-by-moment shell L1 change
is at most 2.24e-6 relative to the original in this selected step,
larger than for the partial lift. The net energy and count
defects are at most 1.39e-14 and 8.22e-15, respectively. The positive exterior fraction
of inverse light products is about 1.2e-11 during the tested kicks.
Refining the common-chart quadrature changes the corrected first-cell
slope by 2.34e-6 fractionally. The maximum corresponding whole-state
L1 change is 8.14e-7 (scalar number), so the larger whole-state
sensitivity itself is not a precision estimate. The first-cell sign is independently corroborated, but the full
common-chart time history and unselected interactions are not tested.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_pair_common_fermion_chart_reference.py --thermal
python -W error warm_pair_common_chart_block_reference.py EOS.dat \
 results/warm-cold-pair-spike-q32-step{10,11,12}-2026-11-26.json.gz \
 results/warm-cold-pair-multishell-q32-9-lifted-2026-11-26.json.gz new-block
~~~

[Matched q32 replay](results/warm-pair-common-chart-q32-replay-2026-11-26.json);
[exact terminal state](results/warm-pair-common-chart-q32-terminal-2026-11-26.json.gz);
[synthetic balance check](results/warm-pair-common-chart-thermal-2026-11-26.json).

A separate 2/4/8-panel test changes only the empty final-fermion
virtual radial grid. Four to eight panels changes the full inverse
moment shapes by at most 1.15e-12 relative to their corresponding
total, but this should not be read as convergence of the tiny
exterior rate itself (1.19e-11 versus 1.20e-11 of inverse product
number). The virtual state is not a formerly untracked input tail.

~~~sh
# Add --fine to the block replay in a fresh directory for the second order.
python -W error warm_pair_chart_exterior_controls_reference.py  results/warm-cold-pair-spike-q32-step12-2026-11-26.json.gz
~~~

[Exact higher-order terminal state](results/warm-pair-common-chart-q32-fine-terminal-2026-11-26.json.gz);
[virtual radial control](results/warm-pair-common-chart-q32-virtual-2026-11-26.json).

## Wider frozen q64 cross-check

At the actual q64 7.754805 MeV native checkpoint, the common
chart and the archived forward projection agree in aggregate across
all four represented/exterior moment arrays to better than 1e-15
in normalized L1 at each matched rule. At the finest tested rule
the global inverse event number and energy match the original
inverse chart to 5.12e-6 and 3.71e-6 relatively. This aggregate
agreement does not test where gain was placed: the first-cell
gain/loss ratio is 5.05, 0.831 and 0.845 on the three original
rules, but about 0.502856 at every tested common rule. The independent
finite tagged integral gives 0.502856094, within 2.9e-7 fractionally
of the finest common value. At that rule the largest normalized
distribution difference is 2.15e-4 for the full signed pair
rate (not for a time-evolved spectrum), down from 3.07e-3 at
the original rule. This comparison freezes the wider state.

~~~sh
python -W error warm_pair_common_fermion_chart_reference.py \
 results/warm-cold-pair-spike-q64-step12-2026-11-26.json.gz
~~~

[Three-rule q64 frozen comparison](results/warm-pair-common-chart-q64-frozen-2026-11-26.json).

## Actual q64 common-chart block

Replay the same q64 late step with the complete common chart active
in all four pair-kick evaluations, not just as a frozen diagnostic.
Use the accepted nine-cell q64 state only as a solver guess. The
independently solved residual falls from 3881.28 to 5.69e-12.
The first-shell occupation slope is −4.71723e-5 per ln(a),
compared with native +1.43554e-4 and nine-cell −4.71724e-5.
The maximum normalized state L1 change from the native trajectory
is 2.17e-6, and the step energy and count defects are about
1.96e-15 and 8.9e-16. This is one accepted wider block at the
original pair quadrature, not the full evolved wider convergence test.

~~~sh
python -W error warm_q64_common_chart_guided_replay_reference.py EOS.dat \
 results/warm-cold-pair-spike-q64-step{11,12}-2026-11-26.json.gz \
 results/warm-cold-pair-multishell-q64-guided-2026-11-26.json.gz common-q64.json
~~~

[Exact q64 output](results/warm-pair-common-chart-q64-terminal-2026-11-26.json.gz)
and [solve diagnostics](results/warm-pair-common-chart-q64-replay-2026-11-26.json).

The actual wider block also completed with the finer (4,10,10,4)
pair rule, changing only that quadrature. The first cold-shell
slope is −4.7172395e-5, within 2.34e-6 fractionally of the coarse
common-chart value. The largest whole-state normalized L1 change
between orders is 8.14e-7, again in the scalar number distribution;
this is not a continuum-precision claim for the whole spectrum.
The guide stalled at scaled residual 8.12e-12 and the original
fallback route accepted the finer result with energy and count
defects 8.39e-15 and 5.33e-15. The guide was not substituted for
a physical state. Reproduce with `--fine` and a separate output path.

[Finer q64 solve diagnostics](results/warm-pair-common-chart-q64-replay-fine-2026-11-26.json)
and [exact finer wider endpoint](results/warm-pair-common-chart-q64-terminal-fine-2026-11-26.json.gz).
