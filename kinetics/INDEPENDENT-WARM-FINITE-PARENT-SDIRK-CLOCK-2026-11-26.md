# Test an L-stable two-stage collision clock without silently changing the physics

On the [common expanding parent chart](INDEPENDENT-WARM-FINITE-PARENT-COMOVING-EXPANSION-2026-11-26.md)
test the stiffly accurate two-stage SDIRK scheme with diagonal
coefficient 1 − 1/sqrt(2). Only the collision subproblem changes;
the centered expansion and the competing pair half-kicks remain
the same. Use the [absolute elastic evaluation](INDEPENDENT-WARM-ABSOLUTE-ELASTIC-EVALUATION-2026-11-26.md)
in this variant to avoid a spurious overflow from relative
thermal coordinates. The established backward-Euler default is
unchanged.

Eliminate the finite parent independently at each implicit stage
with the same conservative local decay weights. The second stage
uses both the daughter-moment and parent affine RK targets made
from the accepted first stage. The affine target itself need not
be an admissible occupation. Both actual stage phases must obey
the Bose and fermion domains, retain the total parent+daughters
invariants, and obey the consistency check on collision entropy.
A failed stage solve is rejected, not silently replaced by
backward Euler or an occupation floor. The scalar stability
function has the correct order conditions and tends to zero in
the stiff limit, but can become negative on the real negative
axis: L-stability alone does not imply positivity.

## Same-phase .68 GeV comparison

The exact .68 n12 phase and the same common parent quadrature,
physical threshold, radiation prescription, interval, pair and
elastic settings are used in all four runs. All final W values
are below the starting value.

| Implicit collision step method | One-block endpoint W / initial daughter energy | Two-block endpoint | Relative endpoint change |
| --- | ---: | ---: | ---: |
| Backward Euler | 2.426473233e-18 | 2.432645582e-18 | +0.2544% |
| SDIRK2 | 2.442816594e-18 | 2.441662125e-18 | −0.04726% |

The measured absolute shift is 5.35 times smaller in this
selected control; that comparison alone is not a measured global
order, a proof of terminal accuracy, or evidence for a new physical
process. Each accepted SDIRK collision-stage invariant defect is
below 7e-16 relative on these runs; the full expansion and
collision ledger is below 3.7e-14 relative.

## Negative .80 GeV stage test

The same SDIRK variant did not find an admissible second-stage
solution on the first block of the tested six-block .80 clock,
although backward Euler accepts the same overall problem. The
affine daughter target left the simple local moment cone in one
cell, and the stage solver stopped near normalized residual 0.43
for three seeds, even after replacing the overflow-prone elastic
form. An affine weight outside the cone is not by itself a proof
that the implicit system has no physical solution. Here it is
a practical warning to retain the verified positive method or
test smaller steps on this .80 phase.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
# Use the exact turn12 and guides080 inputs of the linked clock audits.
d=sdirk
mkdir -p "$d"
run() { python -W error warm_finite_parent_comoving_expansion_reference.py EOS.dat "$@"; }
run "$d/068-abs1" expand-finite-all turn12/step-{8,10}.json --order 2 --common-subdivisions 2 --collision-clock sdirk2
run "$d/068-abs2" expand-finite-all turn12/step-{8,9,10}.json --order 2 --common-subdivisions 2 --collision-clock sdirk2
python -W error warm_coupled_sdirk_audit_reference.py "$d" comoving/068-locked1-expand-finite-all comoving/068-locked2-expand-finite-all
# Separate negative test: this command is expected to reject a stage.
run "$d/080-negative" expand-finite-all guides080/step-{0..6}.json --order 2 --common-subdivisions 2 --collision-clock sdirk2
~~~

[Method comparison and conservation audit](results/warm-finite-parent-sdirk-clock-audit-2026-11-26.json);
[accepted .68 phase records](results/warm-finite-parent-sdirk-hard068-2026-11-26.json.gz).

## .80 GeV short-prefix acceptance and hard-shell rate

A much shorter selected prefix, with Δlog a = 1.875e-8, accepts
both schemes from the same .80 initial phase and common parent
chart. The endpoint W values are 6.390821982e-14 (backward Euler)
and 6.377395203e-14 (SDIRK2). The latter changes by −0.2101%
relative to the former, and both are below the same initial W.
SDIRK2 has stage invariant defects at or below 4.45e-16 relative.
This is a short-prefix result, not an accepted high-order run
over the full original 9e-7 interval.

In a separate frozen ordered-rate diagnostic, the .795–.805 GeV
light shell of that same initial phase has a prompt loss per
particle per Hubble time of 3.97e7, with a gain of 1.81e6.
That is about 5.96 gross-loss turnovers over a 1.5e-7 block
at the initial rate. It helps explain why a longer step can
be difficult, but the turnover is not a linear eigenvalue or
a proof of a positivity timestep threshold.

~~~sh
# First reproduce guides080-12 with the twelve-block recipe.
g=guides080-12
python -W error warm_late_competing_intermediate_guide_reference.py EOS.dat "$g/step-0.json" "$g/step-1.json" half.json
python -W error warm_late_competing_intermediate_guide_reference.py EOS.dat "$g/step-0.json" half.json quarter.json
for m in be sdirk2; do
 run "$d/080-quarter-$m" expand-finite-all "$g/step-0.json" quarter.json --order 2 --common-subdivisions 2 --collision-clock "$m"
done
python -W error warm_prompt_shell_activity_reference.py EOS.dat competing080/tapes/step-{5,6}.json
~~~

[Both accepted .80 short-prefix phases and frozen shell activity](results/warm-finite-parent-sdirk-hard080-short-2026-11-26.json).

## Matched four-block extension for .68 GeV

Start from the same .68 phase and use the identical parent Q grid,
weight fingerprint, energy cut and final scale. The intermediate
native frozen phases only guide the nonlinear roots. Four
accepted blocks produce W = 2.436581812e-18 with backward Euler
and 2.441392797e-18 with SDIRK2. In both clocks W decreases at each
saved step. The absolute two-to-four-block changes are 3.936e-21
and 2.693e-22 of the initial daughter energy, respectively.
For SDIRK the ratio of successive endpoint differences is 4.29,
versus 1.57 for backward Euler. This is consistent with reduced
step sensitivity in this case, not a proof of asymptotic order for
the full split dynamics. The maximum four-block invariant defects
are below 1.66e-13 relative.

~~~sh
# With the existing run function, d=sdirk and the archived turn12 inputs.
python -W error warm_late_competing_intermediate_guide_reference.py EOS.dat turn12/step-8.json turn12/step-9.json "$d/mid89.json"
python -W error warm_late_competing_intermediate_guide_reference.py EOS.dat turn12/step-9.json turn12/step-10.json "$d/mid910.json"
for mode in be sdirk2; do
 dest="$d/068-be4"; [ "$mode" = sdirk2 ] && dest="$d/068-abs4"
 run "$dest" expand-finite-all turn12/step-8.json "$d/mid89.json" turn12/step-9.json "$d/mid910.json" turn12/step-10.json --order 2 --common-subdivisions 2 --collision-clock "$mode"
done
python -W error warm_coupled_sdirk_ladder_reference.py "$d" comoving
~~~

[Matched one/two/four-block audit](results/warm-finite-parent-sdirk-hard068-four-clock-2026-11-26.json);
[accepted four-block phases](results/warm-finite-parent-sdirk-hard068-four-states-2026-11-26.json.gz).
