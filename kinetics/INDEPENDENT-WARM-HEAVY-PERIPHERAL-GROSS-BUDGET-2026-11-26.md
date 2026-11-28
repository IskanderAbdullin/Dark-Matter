# What must be measured beyond the retained daughter support?

The [finite-support comparison](INDEPENDENT-WARM-HEAVY-HIGH-BAND-FINITE-DAUGHTER-SUPPORT-2026-11-26.md)
is conditional. It need not be used as an all-or-nothing
zero-exterior hypothesis. For **number-changing** events
creating a heavy particle above the 1 GeV comoving cut,
at least one dark incoming light/scalar momentum must be
above its retained cap P. If only one is outside, the
energy of the core particle, divided by any high heavy
energy, is at most max(P/P_cut,m_light/M,m_scalar/M).
This follows by comparing both the mass and squared
comoving momentum terms before taking the square root.
Consequently the outside input supplies at least one
minus that ratio of the energy of the high heavy output,
even if both outputs are counted when both are high.

The total comoving energy credited to high heavy states
is therefore at most **gross** outside input to these
events divided by that minimum share. Count an outside
cohort again every time it supplies another such event.
Instantaneous outside energy, or net outward flux, does
not give this gross bound. Loss of high heavy particles
can only improve the deposition comparison; unbounded
heavy elastic transfer across the cut remains excluded.

For q32, the marginal energy share of an outside
participant is at least .6962. With the previous
below-cut and visible allowances already included,
gross peripheral inputs below about 1.41% of the full
comoving budget would suffice for the **partial
direct-light** heavy-energy criterion. This is a target
for a kinetic ledger, not a measured tail bound, and
does not control light decay and elastic gain.

~~~sh
D=results
python -W error warm_heavy_peripheral_gross_budget_reference.py \
 "$D"/warm-heavy-high-band-finite-daughter-support-2026-11-26.json
~~~
