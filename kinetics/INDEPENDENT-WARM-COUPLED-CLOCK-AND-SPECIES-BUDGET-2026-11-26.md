# Check the coupled clock and resolve the species energy exchange

The [comoving finite-parent comparison](INDEPENDENT-WARM-FINITE-PARENT-COMOVING-EXPANSION-2026-11-26.md)
combined expansion, finite prompt and inverse decay, below-pole elastic
scattering, and positive outgoing pair conversion in one local window.
Here compare one and two outer blocks with all those selected reactions
active. Match the starting daughters, the initial parent prescription,
the physical cut, and the exact comoving parent nodes and number weights
between all alternatives. The second block starts from the accepted
first state, including the finite parent and receipts. Nothing re-enters
from outside the chosen support.

| Hard packet | Finite branch W / initial daughter energy, one block | Two blocks | Difference due to clock | Finite parent effect on two-block W |
| --- | ---: | ---: | ---: | ---: |
| .68 GeV | 2.426473233e-18 | 2.432645582e-18 | +0.2544% | +1.20e-7 relative |
| .80 GeV | 2.691195724e-14 | 2.222994659e-14 | −17.3975% | +7.94e-9 relative |

The net decline survives at both resolutions in all four matched
branches. For .80, the 17.4% endpoint shift is direct evidence that
one block is not a converged endpoint prediction. It does not quantify
the remaining error at two blocks.

## Changes split into scalar and light-fermion daughters

Independently integrate number, kinetic energy, total energy and
excess above the fixed physical energy cut from the archived one-block
daughter phases. In the .68 finite fully coupled case, the excess
changes by +1.285e-20 for the scalar and −4.493e-20 for the light
fermion, in units of initial represented daughter energy. The net
decline comes from the light species even as the scalar contribution
rises. In the .80 case both fall, with −3.802e-14 for the scalar
and −1.284e-15 for the light fermion.

The small changes in the *tail* do not imply equally small changes
in the bulk. In the .80 window, adding pair and elastic to finite
prompt shifts the net scalar comoving kinetic energy loss from
−7.20914e-6 to −7.00757e-6 of the initial daughter energy, a 2.80%
reduction in magnitude. In the .68 window the corresponding net
scalar losses are −3.45091e-7 and −3.41683e-7, about a 0.99%
reduction. These are changes of a⁴ times kinetic energy density,
not isolated heat rates or a temperature measurement: expansion
and particle-number exchange are included, and the stated case
comparison does not separately identify elastic and pair causes.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
# Use the same archived turn12 and competing080/tapes inputs as in the linked notes.
d=comoving
run() { python -W error warm_finite_parent_comoving_expansion_reference.py EOS.dat "$@"; }
for mode in expand-stationary-all expand-finite-all; do
 run "$d/068-locked1-$mode" "$mode" turn12/step-{8,10}.json --order 2 --common-subdivisions 2
 run "$d/068-locked2-$mode" "$mode" turn12/step-{8,9,10}.json --order 2 --common-subdivisions 2
 run "$d/080-locked1-$mode" "$mode" competing080/tapes/step-{5,6}.json --order 2 --common-subdivisions 2
 run "$d/080-locked2-$mode" "$mode" competing080/tapes/{step-5,mid-5-to-6,step-6}.json --order 2 --common-subdivisions 2
done
for k in 068 080; do
 s=turn12/step-8.json; [ "$k" = "080" ] && s=competing080/tapes/step-5.json
 for n in 1 2; do
  python -W error warm_finite_parent_comoving_expansion_audit_reference.py "$s" "$d" "$k-locked$n-" expand-stationary-all expand-finite-all >"$d/$k-clock-$n.json"
 done
 python -W error warm_finite_parent_comoving_clock_refinement_reference.py "$d/$k-clock-1.json" "$d/$k-clock-2.json" --all
done
python -W error warm_finite_parent_expansion_species_budget_reference.py turn6/step-4.json "$d"/068-expand-{stationary,finite,stationary-all,finite-all}/step-1.json
python -W error warm_finite_parent_expansion_species_budget_reference.py competing080/tapes/step-5.json "$d"/080-expand-{stationary,finite,stationary-all,finite-all}/step-1.json
~~~

[Audited clocks and species budgets](results/warm-coupled-clock-and-species-budget-2026-11-26.json),
[.68 accepted states](results/warm-coupled-clock-hard068-2026-11-26.json.gz),
[.80 accepted states](results/warm-coupled-clock-hard080-2026-11-26.json.gz).

## Instantaneous pressure and velocity moments

On those same matched one-block endpoints, integrate pressure, the
massive pressure trace and the number-weighted squared speed directly
on the physical momentum grid. Split quadrature at p/m = .1, 1 and 10,
cross-check 32 and 48 nodes, verify 3P + trace = energy, and validate
number and energy independently against the stored shell moments.

Despite the 2.80% change in the .80 scalar kinetic-energy *loss over the
interval*, the scalar endpoint squared speed remains about
0.999999999922 in both models. At fixed endpoint scale, adding pair and
elastic changes the light fermion squared speed from 0.5474095600
to 0.5474095858, a relative shift of 4.71e-8. Its number fraction
below p = m falls by 3.94e-8 in absolute fraction. The .68 example
has still smaller shifts, +7.53e-10 in squared speed in relative units
and −6.58e-10 in fraction below p = m. These are instantaneous
diagnostics of the two selected warm daughter spectra, not a free
streaming integral, a different thermal fit, or a matter-power result.

~~~sh
python -W error warm_expanding_daughter_pressure_reference.py turn6/step-4.json comoving/068-expand-{finite,finite-all}/step-1.json
python -W error warm_expanding_daughter_pressure_reference.py competing080/tapes/step-5.json comoving/080-expand-{finite,finite-all}/step-1.json
~~~

[Species pressure, trace, speed and cold-fraction diagnostics](results/warm-coupled-clock-pressure-moments-2026-11-26.json).

## .80 GeV: actual 3- and 6-block runs on the same comoving chart

The 17.40% one-to-two-block drift persists as a substantial time
resolution issue. Reuse the exact same starting phase, physical
cut and parent quadrature nodes and number weights. Recut the lab
angular geometry on each accepted block. The archived frozen
prompt microstates are used only as solver initial guesses.

| Outer blocks | Finite-parent fully coupled W / initial daughter energy |
| ---: | ---: |
| 1 | 2.691195724e-14 |
| 2 | 2.222994659e-14 |
| 3 | 2.031036969e-14 |
| 6 | 1.814219661e-14 |

Six blocks lower the result by 10.68% from three. The observable
decreases at every saved step in this ladder, while its magnitude
is still visibly sensitive to time resolution. No continuum value
or numerical error bar is inferred. The 1- to 6-block change in
instantaneous full daughter pressure is much smaller in relative
units, −4.52e-6 for the scalar and +3.42e-6 for the light fermion.
Those pressure differences are endpoint diagnostics, not a
cosmological transfer calculation.

The stationary-parent branch also decreases at each saved step.
At six blocks the finite-parent correction to W is +1.01e-08
relative to the matched stationary-parent branch, far smaller
than the measured clock effect.

~~~sh
# Prepare competing080/tapes as above, then regenerate the six exact native guide roots.
t=competing080/tapes
python -W error warm_late_competing_locked_guides_reference.py EOS.dat "$t/step-5.json" "$t/step-6.json" guides080 --steps 6
for mode in expand-stationary-all expand-finite-all; do
 run "$d/080-locked3-$mode" "$mode" guides080/step-{0,2,4,6}.json --order 2 --common-subdivisions 2
 run "$d/080-locked6-$mode" "$mode" guides080/step-{0,1,2,3,4,5,6}.json --order 2 --common-subdivisions 2
done
python -W error warm_expanding_clock_ladder_reference.py "$d" 080 1 2 3 6
parts=(); for n in 1 2 3 6; do parts+=("$d/080-locked$n-expand-finite-all/step-$n.json"); done
python -W error warm_expanding_daughter_pressure_reference.py "$t/step-5.json" "${parts[@]}"
~~~

[Four-step-count ladder and the pressure resolution comparison](results/warm-coupled-hard080-deeper-clock-audit-2026-11-26.json);
[3/6-block accepted daughter phases](results/warm-coupled-hard080-deeper-clock-tapes-2026-11-26.json.gz).

## .80 GeV: matched twelve-block refinement

The same comoving parent chart, prescribed radiation clock, scale
endpoint and initial phase can also be used with twelve outer
blocks. Both branches decrease at every saved step. Their final W
values are 1.6945469363e-14 (stationary parent) and 1.6945469548e-14
(finite parent). Relative to six blocks this is another −6.596%
in either branch, whereas the relative finite-parent change at
twelve blocks is +1.09e-8. The cumulative energy ledger defect
stays below 2.24e-13 relative. This reinforces the distinction
between the robust direction of the decline and an unresolved
continuum amplitude. It does not establish that twelve blocks
are accurate enough.

On the finite-parent branch the one-to-twelve-block pressure
changes remain much smaller than the excess-energy shift:
−5.10e-6 for scalar pressure and +3.86e-6 for light-fermion pressure
in relative units. No matter-power transfer is computed here.

Generate twelve new frozen prompt guide roots on the finer clock.
The guide endpoint need not equal the six-step archive, and the
guides affect only Newton initialization. The prior six-guide
reproduction still requires the original terminal-moment match.

~~~sh
t=competing080/tapes
python -W error warm_late_competing_locked_guides_reference.py EOS.dat "$t/step-5.json" "$t/step-6.json" guides080-12 --steps 12 --refine-clock
for mode in expand-stationary-all expand-finite-all; do
 run "$d/080-locked12-$mode" "$mode" guides080-12/step-{0..12}.json --order 2 --common-subdivisions 2
done
python -W error warm_expanding_clock_ladder_reference.py "$d" 080 1 2 3 6 12
parts=(); for n in 1 2 3 6 12; do parts+=("$d/080-locked$n-expand-finite-all/step-$n.json"); done
python -W error warm_expanding_daughter_pressure_reference.py "$t/step-5.json" "${parts[@]}"
~~~

[Five-clock comparison and pressure diagnostic](results/warm-coupled-hard080-twelve-clock-2026-11-26.json);
[finite-parent accepted states](results/warm-coupled-hard080-twelve-finite-2026-11-26.json.gz);
[stationary-parent accepted states](results/warm-coupled-hard080-twelve-stationary-2026-11-26.json.gz).

## Where the six-to-twelve-step difference is visible

Compare the finite-parent six- and twelve-block W histories at
their six common saved times, not only the endpoint. The coarser
value is larger at each common sample. Their difference is already
1.016% of the starting W at Hdt = 2.50e-7, peaks at about 1.870%
of starting W at Hdt = 8.50e-7, and ends at 1.807% of starting W.
Relative to the *current fine value*, the difference instead
increases from 1.290% to 7.062% across those samples. This locates
substantial step sensitivity before the final interval, but the
recorded difference is cumulative, not a local truncation error.
It does not by itself establish an optimal step schedule.

~~~sh
python -W error warm_clock_path_error_profile_reference.py comoving/080-locked6-expand-finite-all comoving/080-locked12-expand-finite-all
~~~

[Matched-time differences and state audit](results/warm-coupled-hard080-path-clock-audit-2026-11-26.json).
