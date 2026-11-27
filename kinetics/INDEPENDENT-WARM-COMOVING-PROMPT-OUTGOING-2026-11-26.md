# Open a finite outgoing daughter pad during coupled parent evolution

The [earlier frozen pad comparison](INDEPENDENT-WARM-LATE-PROMPT-PAD-OPENING-2026-11-26.md)
allowed outgoing angles while keeping the daughter phase fixed.
Here let particles newly produced outside the original 1.05 GeV
daughter momentum edge evolve and participate in the selected
finite-parent prompt reaction. Include the same local expansion
and parent transport as in the [comoving evolution](INDEPENDENT-WARM-FINITE-PARENT-COMOVING-EXPANSION-2026-11-26.md).
No pair or elastic operator is turned on in this particular
support comparison.

For each tested pad, both branches have the same additional shells,
same exact vacuum occupation there, same parent-grid weights,
same initial parent equilibrium built on the original closed
daughter support, and same endpoint guide for the original shells.
On the control branch the original outer angular boundary stays
closed. The other branch allows on-shell prompt reactions into
the extra shells. A positive predictor in those shells is only
a Newton guess, not an initial occupation floor. All retained
number and energy ledgers include the finite parent.

| Starting packet | Comoving pad edge | Newly active pad energy / initial daughter energy | Newly active pad W / initial daughter energy |
| --- | ---: | ---: | ---: |
| .68 GeV | 1.30 GeV | 3.954e-31 | 1.475e-31 |
| .68 GeV | 1.45 GeV | 4.080e-31 | 1.537e-31 |
| .80 GeV | 1.30 GeV | 3.937e-30 | 1.335e-30 |
| .80 GeV | 1.45 GeV | 4.072e-30 | 1.404e-30 |

All four open cases populate both new species supports in the
chosen pad. The endpoint W integrated over the full resolved
daughter phase is essentially unchanged at this precision:
relative differences from controls are all below 4.3e-14 in the
first four comparisons. Do not infer the sign of an exact
net change at that scale. Independently integrate both new-shell
contributions directly in momentum, check the snapshot and node
fingerprints, and recompute W across the full daughter phase.

For the .80 GeV, 1.30 GeV pad, tighten the nonlinear residual
gate from 8e-12 to 8e-13. Both branches converge below 6.6e-13;
the newly populated pad W is consistent with the default root,
and the relative endpoint comparison is again below 8e-15.
An attempted 3e-13 gate for this open branch stopped at 6.58e-13
when Newton could not reduce its normalized residual further.
No result is claimed at that tighter tolerance.

This is a two-cap check on a conditional finite-volume model. Neither
pad includes arbitrary exterior incoming distributions or guarantees
that all possible subsequent outgoing momenta are represented.
It is not a solution with fully open boundary conditions or a
bound at other times, momenta, or reaction channels.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
# Prepare turn6 and competing080/tapes as in the linked notes.
d=outgoing
mkdir -p "$d"
run() { python -W error warm_finite_parent_comoving_expansion_reference.py EOS.dat "$@"; }
for k in 068 080; do
 paths=(turn6/step-{4,5}.json)
 if [ "$k" = 080 ]; then paths=(competing080/tapes/step-{5,6}.json); fi
 for p in 1.30 1.45; do
  run "$d/$k-$p-closed" expand-finite "${paths[@]}" --order 2 --prompt-pad "$p" --restrict-prompt-pad
  run "$d/$k-$p-open" expand-finite "${paths[@]}" --order 2 --prompt-pad "$p"
  python -W error warm_comoving_outgoing_expansion_audit_reference.py "${paths[0]}" "$d/$k-$p-closed" "$d/$k-$p-open"
 done
done
p=1.30
for mode in closed open; do
 opts=(); [ "$mode" = closed ] && opts=(--restrict-prompt-pad)
 run "$d/tight-$mode" expand-finite competing080/tapes/step-{5,6}.json --order 2 --prompt-pad "$p" --root-tol 8e-13 "${opts[@]}"
done
python -W error warm_comoving_outgoing_expansion_audit_reference.py competing080/tapes/step-5.json "$d/tight-closed" "$d/tight-open"
~~~

[Independently audited pad receipts](results/warm-comoving-prompt-outgoing-audit-2026-11-26.json);
[.68 accepted daughter phases](results/warm-comoving-prompt-outgoing-hard068-2026-11-26.json.gz);
[.80 accepted daughter phases](results/warm-comoving-prompt-outgoing-hard080-2026-11-26.json.gz).
