# Separate pair conversion from elastic energy transfer

The [earlier coupled comparisons](INDEPENDENT-WARM-COUPLED-CLOCK-AND-SPECIES-BUDGET-2026-11-26.md)
turned on pair and elastic channels together. Here keep finite
parent transport, prompt decay and inverse decay, and the same
radiation clock in all cases. For each archived phase and
discretization, compare four equal-endpoint versions: prompt
only, prompt plus pair, prompt plus elastic, and all three.
Within each comparison the parent quadrature, physical cut,
EOS and accepted block times match exactly. The one- and
two-block comparisons use their own previously validated
charts; the .68 tests also have distinct archived phases.
Their agreement is a repeat across these controls, not a
pure isolation of the time discretization.

Positive numbers in the table mean an increase in final scalar
a⁴ times kinetic energy density relative to the prompt baseline,
in units of initial total daughter energy density. The prompt
baseline still loses net scalar kinetic energy in each case.

| Packet | Blocks | Pair increment | Elastic increment | Both |
| --- | ---: | ---: | ---: | ---: |
| .68 GeV | 1 | +3.832e-9 | −4.240e-10 | +3.408e-9 |
| .68 GeV | 2 | +3.826e-9 | −4.242e-10 | +3.402e-9 |
| .80 GeV | 1 | +2.266e-7 | −2.506e-8 | +2.016e-7 |
| .80 GeV | 2 | +2.250e-7 | −2.508e-8 | +1.999e-7 |

The pair contribution partly offsets the prompt scalar loss.
Elastic scattering alone moves energy in the opposite direction,
increasing light-particle kinetic energy at the expense of the
scalar in these controls. The separate scalar increments
numerically almost add: the residual of the factorial contrast
is below 2.7e-13 in these units. This is a discrete comparison
of final states, not a statement that the channels are independent.

About 60–61% of the pair-only positive scalar kinetic increment is
the increased comoving scalar count evaluated at the prompt-only
mean. The rest is the change in comoving mean kinetic energy
evaluated at the new count. This exact algebraic decomposition
separates count from mean-energy change without assigning a
temperature to a nonthermal phase. The full data also record
kinetic/total energy, number, physical excess W and daughter
entropy separately. Daughter entropy alone is not a collision
entropy balance for the entire system.

At .80 in one block, pair lowers the total W relative to the
prompt comparison by 6.77e-20 in units of initial energy and
elastic lowers it by 1.90e-21. Thus the small sensitivity of W
does not imply the kinetic ledger is insensitive. No observable
large-scale transfer function is inferred here.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
# Prepare turn6, turn12 and competing080/tapes as in the linked notes.
d=factorial; b=comoving
mkdir -p "$d"
run() { python -W error warm_finite_parent_comoving_expansion_reference.py EOS.dat "$@"; }
a() { python -W error warm_coupled_channel_factorial_audit_reference.py "$@"; }
for mode in expand-finite-pair expand-finite-elastic; do
 run "$d/068-n1-$mode" "$mode" turn6/step-{4,5}.json --order 3 --common-subdivisions 1
 run "$d/068-n2-$mode" "$mode" turn12/step-{8,9,10}.json --order 2 --common-subdivisions 2
 run "$d/080-n1-$mode" "$mode" competing080/tapes/step-{5,6}.json --order 3 --common-subdivisions 1
 run "$d/080-n2-$mode" "$mode" competing080/tapes/{step-5,mid-5-to-6,step-6}.json --order 2 --common-subdivisions 2
done
run "$d/080-n2-expand-finite" expand-finite competing080/tapes/{step-5,mid-5-to-6,step-6}.json --order 2 --common-subdivisions 2
a turn6/step-4.json "$b/068-expand-finite" "$d/068-n1-expand-finite-pair" "$d/068-n1-expand-finite-elastic" "$b/068-expand-finite-all"
a turn12/step-8.json "$b/068-locked2-expand-finite" "$d/068-n2-expand-finite-pair" "$d/068-n2-expand-finite-elastic" "$b/068-locked2-expand-finite-all"
a competing080/tapes/step-5.json "$b/080-expand-finite" "$d/080-n1-expand-finite-pair" "$d/080-n1-expand-finite-elastic" "$b/080-expand-finite-all"
a competing080/tapes/step-5.json "$d/080-n2-expand-finite" "$d/080-n2-expand-finite-pair" "$d/080-n2-expand-finite-elastic" "$b/080-locked2-expand-finite-all"
~~~

[Four controlled channel comparisons](results/warm-coupled-channel-factorial-audit-2026-11-26.json);
[.68 accepted cases](results/warm-coupled-channel-factorial-hard068-2026-11-26.json.gz);
[.80 accepted cases](results/warm-coupled-channel-factorial-hard080-2026-11-26.json.gz).

## Independent frozen pair rate below and above the cut

Project the positive forward and reverse pair contributions
separately onto the exact cut-aligned archived shell edges.
For the .80 phase at Hdt = 1e-7, scalar net kinetic gain below
the cut is about +0.25695 in units of total daughter energy
per Hubble time. Above the cut, the net scalar kinetic current
is negative, about −1.8e-12 to −1.9e-12 in the same units in two
quadrature orders. The .68 phase similarly has positive
below-cut scalar gain, about +0.25694. This independent static
fold locates the scalar energy source on the bulk side of the
cut; it is not an integral of the evolving energy ledger.

~~~sh
python -W error warm_frozen_pair_species_regions_reference.py EOS.dat turn6/step-4.json competing080/tapes/step-{5,6}.json
~~~

[Frozen directional pair rates, split by species and energy region](results/warm-frozen-pair-species-regions-2026-11-26.json.gz).
