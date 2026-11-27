# Open exterior receipts from late below-pole elastic scattering

The matched late [.68](INDEPENDENT-WARM-LATE-COMPETING-FLUX-2026-11-26.md)
and [.80 GeV](INDEPENDENT-WARM-LATE-HARD080-COMPETING-FLUX-2026-11-26.md)
tests evolved closed elastic scattering. Here hold each set of three
saved daughter phases fixed, keep zero incoming occupation above
the original 1.05 GeV edge, and open all possible high outgoing
angles in the specified below-pole channel (s < 0.82 M²). No exterior
particles are returned to the reacting phase. The largest possible
incident lab energy from these two supports is about 2.1002 GeV, so
a pad extending to that momentum cannot remove a possible outgoing
event at the high edge. The original low boundary is unchanged.

On the expanded chart the original closed current is measured with
exactly the same quadrature samples, by excluding events leaving the
old domain. As a separate control it is also recomputed on the old
chart with its own breakpoints. For every populated cross-support
event, exactly one outgoing daughter is counted in the exterior, as
expected from the total energy ceiling. The number and energy of
those receipts reproduce W outside independently via E − C N,
where C = 0.7148213914435242 GeV. We sum the very small crossing
terms directly rather than subtracting total collision currents.

In six tested phases and two orders, the largest sampled exterior
receipt is 4.67e-30 in units of initial dark energy per Hubble time.
The maximum of the receipt and the absolute *net W change* in
cross-support events is 1.67e-19 of the absolute prompt current.
These are observed maxima, not mathematical upper bounds. The
exterior receipt itself is still quadrature-sensitive: for example,
the first .68 value moves from 4.84e-31 to 2.15e-30 when increasing
the order. That does not change the scale comparison on these
samples. This is not a new time integration or a constraint on
incoming omitted populations or higher-pole reactions.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
# Prepare turn6 and competing080/tapes as in the linked notes.
for mode in 0 1; do
  python -W error warm_late_open_elastic_receipt_reference.py EOS.dat "open-$mode.json"     turn6/step-{0,4,5}.json competing080/tapes/step-{0,5,6}.json --order-mode "$mode"
done
python -W error warm_late_open_elastic_receipt_audit_reference.py open-0.json open-1.json
~~~

[Ledger audit](results/warm-late-open-elastic-receipt-audit-2026-11-26.json);
[two raw quadrature tapes](results/warm-late-open-elastic-receipt-tapes-2026-11-26.json.gz).
