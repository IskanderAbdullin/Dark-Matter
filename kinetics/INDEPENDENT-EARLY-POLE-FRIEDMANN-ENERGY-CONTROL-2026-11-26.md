# Live parent and separate conserved visible radiation in Friedmann H

An independent follow-up applies the
[energy-conserving collision stencil](INDEPENDENT-EARLY-POLE-EXACT-ENERGY-PARENT-STENCIL-2026-11-26.md)
to a separated visible radiation bath and a dark parent plus two dark
daughters. Radiation has constant visible degrees of freedom
g* = 106.75 and no heat exchange with the dark sector, so its
temperature remains proportional to 1/a. The active dark species all
contribute to H through their energy density. The clock normalization
is defined against visible radiation alone, avoiding double counting
these explicit populations in g*.

The variational Hartree energy supplies its bare-mass trace, not a
replacement trace proportional to the self-induced pole. The massive
parent contributes its integrated shell trace. Integrating the trace
independently across BDF steps checks the change in total comoving
energy, not only the instantaneous collision identity.

In the refined case, tracked dark energy is 3.262% of visible
radiation; H is 1.618% greater than the reference. The yield by x=0.05
changes from 0.000350663 (radiation-only H) to 0.000346891 (full H),
a decrease of 1.08%. On the refined energy-evolving run the total
continuity residual is less than 10⁻¹¹ of the final daughter energy;
the coarse residual is 1.23 × 10⁻¹¹. These figures are conditional on
this initial thermal parent, effective dark interactions, constant
g*, and interval. They do not include reheating, chemical feeding,
or a primordial/late-time cosmological fit.

Combining isolated energy-conserving parent mixing with the
conservative collision stencil gives yields 0.000352298 and
0.000406585 at rates 3 and 100 per reference expansion time with
radiation-only H, compared with 0.000350663 without mixing. This
reproduces the partial recovery of scalar yield with fast mixing.
With the same rate 100 and full H the yield is 0.000402572.
The slow-rate change was comparable to discretization shifts in the
older controls and does not establish a definite sign. The mixing
rate is phenomenological, not a scattering-rate calculation.

~~~sh
python -W error early_pole_friedmann_energy_control_reference.py
python -W error early_pole_friedmann_energy_control_reference.py --scalar 65 --fermion 90 --parent 72
python -W error early_pole_friedmann_energy_control_reference.py --scalar 65 --fermion 90 --parent 72 --radiation-only-H
python -W error early_pole_elastic_energy_combined_reference.py --rate 3
python -W error early_pole_elastic_energy_combined_reference.py --rate 100
python -W error early_pole_elastic_energy_combined_reference.py --rate 100 --full-H
~~~

[Raw comparisons](results/early-pole-friedmann-energy-control-2026-11-26.json).
