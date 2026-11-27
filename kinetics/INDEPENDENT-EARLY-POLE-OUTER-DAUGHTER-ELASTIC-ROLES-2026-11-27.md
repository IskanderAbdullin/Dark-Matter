# First-window outside-daughter preparation: matched elastic controls

This is a conditional calculation at the same (3.33333\times10^{-7})
joined time window and 16 macro steps, not a channel decomposition of a physical
observable. For each of four elastic controls, compare the integrated cold
source for the native hard state with the source when only the daughter
moments outside the initial cold aperture are replaced. The comparison checks
exact agreement of hard grids, time schedule and initial cold data and requires
all eight evolutions to pass the existing energy and entropy gates.

| Elastic control | Difference in mean cold source (outside minus native) |
|---|---:|
| None | -0.000000471040 |
| Unlike only | -0.000752345519 |
| Same only | -0.000140525456 |
| Both | -0.000741836280 |

The unlike-only control reproduces the effect to within 1.42% of the two-channel
magnitude; the same-only effect is 18.94% of that magnitude. Effects cannot be
added: the four-arm interaction is +0.000150563655, opposing the single-channel
effects. Full screening feedback is enabled where the unlike packet is on; it is
not a separately toggled factor here. Thus the unlike control identifies
sufficiency in this local model test, not an additive share or a universal
ranking.

An independent initial-state test projects the held elastic moment derivative
onto the four initially cold-active hard cells (two per charged species) and
recomputes the dressed cold source in both directions. At probe amplitudes
(3\times10^{-15}) and (10^{-15}), the total two-channel slope difference
between the prepared and native cases is about (-9.3\times10^3) and
(-9.5\times10^3), respectively. Restricting the derivative to the four active
cells reproduces the slope difference within about 0.2% at these probes.
The directly projected outside contribution is not numerically resolved at
arbitrary precision. This test localizes where the instantaneous derivative
enters the source in the held graph; it does not assign unique microscopic
paths or predict finite-time transport from the initial derivative alone.

The finite receipt is reproduced with the native and outside raw split outputs
in the order none, unlike-only, same-only, both, using
`early_pole_outer_daughter_elastic_factorial_reference.py`.
The direction receipt uses the archived paired bundle with
`early_pole_outer_daughter_initial_flux_locality_reference.py`.
