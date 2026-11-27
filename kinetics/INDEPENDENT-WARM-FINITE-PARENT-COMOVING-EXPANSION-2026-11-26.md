# Compare finite heavy-parent feedback with and without local expansion

Use the [common comoving parent chart](INDEPENDENT-WARM-COMOVING-PARENT-CHART-2026-11-26.md)
to keep the evolving parent attached to its exact momentum nodes.
For each test, start both stationary and finite-parent branches from
the same archived daughter phase. The finite branch starts with the
resolved instantaneous parent balance. At the first half-redshift,
the finite occupation stays fixed at each comoving Q. Then evolve the
parent and daughters together at the midpoint collision chart with
one condensed implicit solve; redshift the accepted phase again.
The stationary comparison uses the same trial parent angle chart but
sets its occupation algebraically at each collision evaluation.

Use a prescribed visible entropy clock and the common passive
initial-parent energy contribution to H. This isolates the effect of
parent time evolution in the reaction step, rather than solving for
a different expansion history in each branch. The heavy parent is
included in the finite branch number and energy ledgers, and its
massive trace work is included across the two drifts. The stationary
branch closes its separate daughter-only ledger. Recompute W with
the fixed physical cut inside the final daughter bins, never with
unchanged cut-aligned weights after redshift.

| Hard packet | Δlog a in expanding branch | Initial W / initial daughter energy | Stationary expanding endpoint | Finite expanding endpoint | Relative endpoint shift |
| --- | ---: | ---: | ---: | ---: | ---: |
| .68 GeV | 1.5e-8 | 2.427745338e-18 | 2.395662461e-18 | 2.395662725e-18 | +1.10e-7 |
| .80 GeV | 9.0e-7 | 6.621907064e-14 | 2.691217331e-14 | 2.691217347e-14 | +6.27e-9 |

Both cases decline on these local clocks. Holding the scale exactly
constant but applying the same collision clock on the same common
chart gives relative finite-parent shifts of +1.10326e-7 and
+6.26924e-9. Changing from the frozen to expanding block shifts
the absolute endpoint more than it shifts the finite correction.
On the independently archived .68 12-step phase, the matched
one- and two-block comparison with *the same comoving weights*
gives relative finite-parent shifts of +1.08e-7 and +1.20e-7,
while changing the number of blocks moves the endpoint about
0.254%. This controls a discrete local conclusion, not a bound on
the exact continuum crossing time.

Allowing the selected open-outgoing pair reaction and closed subpole
elastic scattering together with finite parent transport leaves
the endpoint declines in both cases. Relative to the matched
stationary parent with those same two channels, the finite shifts
are still +1.10e-7 and +6.27e-9. Compared with the finite prompt-only
runs, the extra reactions themselves change the endpoints by
−1.37e-8 and −2.59e-6. Thus the sign result is not an artifact of
running these local corrections only in isolation. Including the
pair receipts and the known species exchange, the combined number
and energy ledgers close to 4.3e-14 relative in the sampled blocks.
The outgoing pair receipts do not feed back into the occupations.

For the .68 chart the radial order-4 endpoint moves by 1.35e-9
relative to order 3 in both branches. Parent-local angular panels
can be generated one at a time instead of being kept in memory;
stored and streamed values and directional derivatives agree
bit-for-bit on the exact archived-phase control.

The model omits visible sources, unknown returning exterior phases
and scalar number-changing processes. This is not a late-time
cosmological integration, a bound outside the chosen interval, or
a measurement of the dark-matter abundance.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
# Prepare turn6, turn12 and competing080/tapes using the cited notes.
d=comoving
mkdir -p "$d"
run() { python -W error warm_finite_parent_comoving_expansion_reference.py EOS.dat "$@"; }
audit() { python -W error warm_finite_parent_comoving_expansion_audit_reference.py "$1" "$d" "${@:2}"; }
modes=(frozen-stationary frozen-finite expand-stationary expand-finite expand-stationary-all expand-finite-all)
for mode in "${modes[@]}"; do
  run "$d/068-$mode" "$mode" turn6/step-{4,5}.json --order 3
  run "$d/080-$mode" "$mode" competing080/tapes/step-{5,6}.json --order 3
done
audit turn6/step-4.json 068- "${modes[@]}"
audit competing080/tapes/step-5.json 080- "${modes[@]}"
for mode in expand-stationary expand-finite; do
  run "$d/068-locked1-$mode" "$mode" turn12/step-{8,10}.json --order 2 --common-subdivisions 2
  run "$d/068-locked2-$mode" "$mode" turn12/step-{8,9,10}.json --order 2 --common-subdivisions 2
  run "$d/068-order4-$mode" "$mode" turn6/step-{4,5}.json --order 4
done
for n in 1 2; do audit turn12/step-8.json "068-locked$n-" expand-stationary expand-finite; done
audit turn6/step-4.json 068-order4- expand-stationary expand-finite
python -W error warm_comoving_streamed_panel_audit_reference.py turn6/step-4.json
~~~

[Verified matched response and conservation](results/warm-finite-parent-comoving-expansion-audit-2026-11-26.json);
[.68 accepted single-block states](results/warm-finite-parent-comoving-expansion-hard068-2026-11-26.json.gz);
[.80 accepted single-block states](results/warm-finite-parent-comoving-expansion-hard080-2026-11-26.json.gz);
[.68 clock and radial controls](results/warm-finite-parent-comoving-expansion-controls-2026-11-26.json.gz).
