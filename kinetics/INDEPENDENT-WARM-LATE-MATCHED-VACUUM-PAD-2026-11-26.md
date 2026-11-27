# Compare late reactivation with the same prompt clock but a vacuum pad

For each of the .68 and .80 GeV test realizations, start
a second run using the identical finite daughter mesh,
fixed H and nominal three-microstep reporting clock.
Keep the full original q64 retained occupation, but do
not activate the saved exterior cohorts in this control.
Their count is carried separately. Both solves use the
same prescribed H from the active realization, rather
than quietly changing the duration. The outside
momentum cells begin at physical vacuum in the control.
A separate locked replay of active .80 intervals 1, 4,
and 6 matches their archived moments to 4.6e-16 or better,
so the comparison does not rely on adaptive path
equivalence at those tested intervals.

The convex excess above the cut is much larger in
the reactivated tests than in their vacuum controls.
At Hdt=1e-6 it is 2.70e-19 of the initial energy for
the .68 case and 2.03e-14 for the .80 case. The matching
vacuum values are 2.86e-43 in both runs. Despite that
large *relative* high-tail change, the largest normalized
shell L1 difference in number, full energy or pressure
at any saved time is below 7.84e-11 in either run.
The restored finite receipt energy is itself very small.
The light mean squared speed shifts by less than
1.51e-11 relative over these outputs. The within-mesh
fourth momentum moment is more tail-sensitive (the
largest scalar shift is 5.65e-8), but this is not an
estimate of a physical transfer function.

## Check the whole reconstructed phase, not just shell moments

At the terminal matching points integrate the absolute
difference of the two analytic occupancies with the
radial measure in each cell. The sign can change only
where the two affine entropy profiles cross. Split
there and compare three radial orders; use stable
exponential differences instead of subtracting two
nearly equal bulk occupancies.

Independently enclose the positive difference on
short subpanels. On each, the minimum entropy sets
the largest admissible Bose/Fermi derivative, and
both endpoints bound the absolute affine entropy
difference. Multiply by the exact radial volume and
sum this conservative upper. Doubling the resulting
unnormalized allowance also accounts for separate
number normalization. For the whole normalized
terminal phase, the L1 upper is below 1.0e-10 for
each species in both tests, up to numerical evaluation
of the positive formula. Any single-particle
observable bounded by one in absolute value then has
an expectation shift at most that L1, for these
particular reconstructed phases. A sampled instantaneous
isotropic velocity characteristic is only an illustration,
not a coupled cosmological Boltzmann hierarchy.

This paired control isolates a finite known exterior
choice in the selected prompt-only test. It cannot
constrain an unrecorded hot distribution, continuous
external sourcing, or the whole early-to-late history.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/warm-corrected-anchor-signed-floor-2026-11-26.json
b=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
s=results/warm-late-cut-aligned-hard068-test-step0-2026-11-26.json.gz
t=results/warm-late-aligned-hard080-test-step0-2026-11-26.json.gz
python -W error warm_late_prompt_matched_vacuum_pad_reference.py EOS.dat "$a" "$b" "$s" vacuum068
python -W error warm_late_prompt_matched_vacuum_pad_reference.py EOS.dat "$a" "$b" "$t" vacuum080
# After reproducing the corresponding three-microstep active runs:
python -W error warm_late_prompt_vacuum_difference_reference.py "$b" active068 vacuum068
python -W error warm_late_prompt_vacuum_difference_reference.py "$b" active080 vacuum080
python -W error warm_late_active_clock_provenance_reference.py EOS.dat active080 1 4 6
python -W error warm_late_prompt_bounded_observable_reference.py active080/step-6.json vacuum080/step-6.json
~~~

[.68 matched series](results/warm-late-vacuum-pad-hard068-compare-2026-11-26.json.gz)
and [.80 matched series](results/warm-late-vacuum-pad-hard080-compare-2026-11-26.json.gz).
Terminal vacuum states: [.68](results/warm-late-vacuum-pad-hard068-terminal-2026-11-26.json.gz)
and [.80](results/warm-late-vacuum-pad-hard080-terminal-2026-11-26.json.gz).
Phase measures: [.68](results/warm-late-vacuum-observable-hard068-2026-11-26.json)
and [.80](results/warm-late-vacuum-observable-hard080-2026-11-26.json).
[Active clock replay](results/warm-late-vacuum-pad-active-clock-control-2026-11-26.json).
