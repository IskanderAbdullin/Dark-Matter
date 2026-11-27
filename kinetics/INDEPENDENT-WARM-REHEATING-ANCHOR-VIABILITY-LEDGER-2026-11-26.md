# What an altered 80 keV history would have to change

This is a feasibility ledger for the selected
[late stable floor](INDEPENDENT-WARM-REHEATING-CAP-FREE-ANCHOR-THRESHOLDS-2026-11-26.md)
and [scalar ceiling](INDEPENDENT-WARM-REHEATING-POST-ANCHOR-SCALAR-CEILING-2026-11-26.md).
Neither inequality alone supplies a viable cosmology.
Keep the stable and scalar requirements separate, and
treat the scenarios below as conditional moment edits
rather than rate calculations.

## Stable depletion before the anchor

Let a different early history deliver only a fraction
$q$ of the nominal 80 keV stable count. Set the altered
late energy budget to
$\rho_{\rm cap}(q)=\rho_{\rm cap}(1)-\eta(1-q)m_\chi n_0$,
where $\eta$ is the fraction of the removed stable rest
energy exported from the dark sector. Thus $\eta=0$
retains all energy for possible dark redistribution,
and $\eta=1$ removes that rest energy from the budget.
This is a specified family of inputs; it says nothing
about the radiation or other fate of the exported energy.
The late comparison still allows the stated lepton input.

Use the endpoint-slab passage lower bound for every
$q$ in a small interval. The positive rate coefficients
are increased to their values at the largest $q$, and
the budget divided by the stable rest energy to its
value at the smallest $q$. The proposed endpoint
$X=\Omega_{\rm goal}/(\Omega_{\rm freeze}q)$ is largest
at the smallest $q$. This avoids treating a spot check
of an individual edited anchor as a result for all
initial counts. Repeat with 256 and 1024 $q$-intervals.

| Fraction of removed rest energy exported | Earlier stable count reduction needed to escape the comparison for \(\Omega_\chi h^2\le .120\) |
|---:|---:|
| 0 | About 56.4% |
| 1/2 | About 80.1% |
| 1 | About 99.4% |

The .060 stable target yields nearly the same
comparison thresholds. These are necessary losses
in the chosen budget family; retaining less than
the indicated initial count may only make the
comparison inconclusive, rather than achieve the
stable target. Export generally makes the survival
obstruction harder to escape because it leaves less
energy available for later stimulated reactions.

## The independent scalar count and entropy problem

The original 80 keV scalar ceiling grants every
remaining weighted particle count to scalars and
includes the allowed charged-lepton input and dark
heavy regeneration. It permits a scalar yield at
most \(6.5936\times10^{-4}\), against
\(2.5285\times10^{-3}\) for \(\Omega_\phi h^2=.060\).
The missing yield is at least \(1.8691\times10^{-3}\)
in that comparison. Simply relabeling stable
particles as scalars is already included in this
upper bound. This scalar ledger uses the original
count and channel assumptions, not an arbitrary
new pre-anchor production history.

As a separate hypothetical change, suppose that
all selected reactions finish before a universal
visible entropy dilution, without new productions
or dark number changes during the dilution. Reducing
the stable lower bound to .120 requires dilution
by at least 1215.5; for .060 it requires 2431.1.
The already insufficient scalars are also diluted:
their upper bounds become \(1.29\times10^{-5}\)
and \(6.44\times10^{-6}\) in \(\Omega h^2\), respectively.
This is not a bound on entropy injection while
the reactions are active or on a simultaneous
new source of scalar particles.

~~~sh
python -W error warm_prompt_reheating_anchor_viability_ledger_reference.py \
 EOS.late.dat \
 results/warm-prompt-reheating-cap-free-anchor-thresholds-2026-11-26.json.gz \
 results/warm-prompt-reheating-post-anchor-scalar-ceiling-2026-11-26.json \
 results/warm-prompt-reheating-cold-lepton-energy-allowance-2026-11-26.json
~~~

[Early depletion and late dilution ledgers](results/warm-prompt-reheating-anchor-viability-ledger-2026-11-26.json).
