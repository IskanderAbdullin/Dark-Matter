# Closing the heat budget with a finite visible radiation bath

This test replaces the [fixed-temperature momentum restoration](INDEPENDENT-EARLY-POLE-LIVE-PARENT-ENERGY-LEDGER-2026-11-26.md)
by a finite, evolving radiation bath, on top of the
[conservative massive-daughter stencil](INDEPENDENT-EARLY-POLE-MASSIVE-DAUGHTER-ENERGY-2026-11-26.md).
Parent momentum relaxation targets a discrete FD shape at the
current bath temperature and at the actual parent number. Whatever
energy it removes from the parent is added to the radiation.
Friedmann H uses the instantaneous energy of all components.
The rate 100 per radiation-only reference expansion time is
illustrative, not a computed scattering rate.

At mass ratio 0.06, with initial visible degrees of freedom 106.75,
the refined bath gains 0.04260 in energy units of the reference
temperature to the fourth power and warms by 0.03031% relative to
free redshifting. This is large relative to the daughter energy but
small relative to the radiation reservoir. Yield relative to reference
entropy is 0.000440067 with finite radiation and 0.000440213 with
fixed bath temperature, a change of −0.033%. Relative to the actual
heated visible entropy, the finite-bath yield is 0.000439667.
The coarse pair shows the same sign and a similar fractional change.

Number balances and a separate continuity integral pass; the refined
energy continuity mismatch is about 10⁻¹¹ of the final daughter
energy. The thermostat comparator logs exported heat without feeding
it back into its temperature or H. This is not a constraint on
heating under an arbitrary different abundance, scattering rate,
or bath heat capacity.

~~~sh
python -W error early_pole_finite_radiation_heat_reference.py
python -W error early_pole_finite_radiation_heat_reference.py --scalar 65 --fermion 90 --parent 72
python -W error early_pole_finite_radiation_heat_reference.py --scalar 65 --fermion 90 --parent 72 --thermostat
~~~

[Raw comparisons](results/early-pole-finite-radiation-heat-2026-11-26.json).
