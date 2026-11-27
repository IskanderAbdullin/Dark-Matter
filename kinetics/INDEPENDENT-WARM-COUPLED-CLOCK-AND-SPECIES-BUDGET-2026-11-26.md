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
