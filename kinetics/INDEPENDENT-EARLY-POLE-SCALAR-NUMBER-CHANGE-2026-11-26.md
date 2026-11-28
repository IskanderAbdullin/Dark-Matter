# Number-changing scalar reactions are distinct from infrared redistribution

For an unbroken real scalar with interaction −λφ⁴/4!, the cold
vacuum threshold 4→2 process has ten tree diagrams. At four incoming
particles at rest, the two equal-energy products have E=2m and
p=√3m. Four internal denominators are 8m² and six are −4m²,
so the signed sum is −1/m² before multiplying by λ². Integrating
two-body phase space and including 4!2! once gives the threshold
**event** coefficient √3 λ⁴/(12288 π m⁸), before multiplying by
the fourth power of number density. Net number loss is twice the
forward-minus-reverse event rate. This agrees with
[Arcadi, Lebedev, Pokorski and Toma](https://arxiv.org/abs/1906.07659),
Eq. (55) and App. C, in that convention.

For a distinct **formal** frozen test, insert the Hartree pole
as a stable mass and the bare threshold amplitude. Restrict the
reaction family to four cold incoming particles with q≤10⁻⁵
and the corresponding two outgoing modes. Random isotropic
cold samples reconstruct the exact center-of-mass two-body
kinematics, all ten propagators and the initial energy factors.
Exact kinematic rates are within 0.01% of the rest threshold
before final occupancy factors. For .06 at x=0.05, the allowed
outgoing q lies near 0.00148–0.00152: outside the very narrow
cold peak even though q is small relative to visible temperature.
The outgoing particles are not nonrelativistic relative to the
Hartree mass.

Using a reference parent mass of 86 GeV, the forward threshold
number-loss clock with empty outgoing states is about 19–20 Href
at x=0.05, and 71–72 Href for the .06 case at x=0.10.
Evaluating both forward and inverse Bose factors in the saved
nonthermal occupations changes the .06 third-grid formal net
number loss to about 70 and 420 Href, respectively. On the
previous grid these become 78 and 315 Href. The sampled hard
shell occupancies, unlike the threshold matrix element, are
important and mesh-sensitive. This is not a rate bound:
higher-order vertices, coherent correlations, spectral widths,
or other families of 4↔2 events can change the result. The
bare threshold treatment is uncontrolled at the snapshot
occupancies, and the previously tested elastic screening factor
cannot be assigned to this different reaction.

~~~sh
python -W error early_pole_scalar_number_change_threshold_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

[Threshold, jump bounds and raw snapshot tests](results/early-pole-scalar-number-changing-threshold-2026-11-26.json).
