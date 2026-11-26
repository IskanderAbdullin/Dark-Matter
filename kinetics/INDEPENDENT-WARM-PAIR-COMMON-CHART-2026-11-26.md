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
(E_{L,max}=2E_{phi,max}-m_L), limits additional virtual fermion
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
