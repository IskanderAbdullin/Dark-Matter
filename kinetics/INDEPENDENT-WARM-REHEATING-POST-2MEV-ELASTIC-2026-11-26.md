# Track chemical depletion and momentum kicks separately below 2 MeV

Continue the [elastic 2 MeV handoff](INDEPENDENT-WARM-REHEATING-ELASTIC-FEEDBACK-2026-11-26.md)
with the same two quantum thermodynamic states and
MB event and energy coefficients. Tabulate the
below-pole heat on the resulting path and re-evaluate
it once. Factor exp(−2m_light/T) analytically out
of the invariant table: the light forward rate then
uses the square of its phase factor at rest, even
when naive quadrature would underflow. Keep the
inverse channel. This is a selected two-temperature
diagnostic, not a kinetic relic calculation.

Both calibrated paths reduce the light yield from
3.19786e−5 at 2 MeV to about 3.0364e−5 at visible
80 keV. The integrated count difference between the
8- and 13-state calibrations is about 2.5e−5 relative.
At the latter sample, the fractional light number-loss
rate is about 3.38e−5 per Hubble time. The level of
momentum redistribution is a different question.
In the finer held calibration, the squared-kick
moment divided by the total light momentum second
moment and H drops below 1 between about 166 and
131 keV. At 80 keV it is about .436 for light and
1.70 for scalars. These are positive collision moments,
not eigenvalues or a demonstration that the thermal
shapes remain valid there.

## A conditional stable remaining-count floor with arbitrary elastic redistribution

An earlier [loss-only ballistic argument](INDEPENDENT-WARM-REHEATING-TWO-TEMPERATURE-LATE-2026-11-26.md#two-independent-conditional-endpoint-checks)
cannot be used if subsequent elastic scatterings
redistribute momentum. Instead use an all-invariant
ceiling on annihilation. The spin-trace Cauchy
comparison bounds the vacuum χχ̄→φφ rate by
y⁴/(8π M²); the lab Møller factor does not increase
the CM limit. Dropping final scalar stimulation in
the vacuum trace is reversed by multiplying this
ceiling by (1+F)² if the occupation of *every*
outgoing scalar mode is at most F at every later
time. The inequality is on the actual density of
both light charges and requires no momentum
restriction or shape assumption. Light inverse
production, if present, adds to that count.

With n_light both charges, nχ nχ̄≤n_light²/4.
If b is the scale factor divided by its anchor
value, the photon expansion floor H≥Hγ/b² and
the factor b⁻³ in the count equation give a
comoving exposure ≤1/Hγ. The resulting retained
fraction is ≥1/[1+C(1+F)²], with
C=n_light(σv)_vac,upper/(2Hγ).
The assumed photon temperature does not fall
below its collisionless b⁻¹ continuation.
Use conserved visible entropy for the rest-density
normalization.

An optional comparison adds all three charged-lepton
final channels with the specified finite vacuum
vector width. Bound s/[(s−V²)²+(VΓ)²] by its
global maximum, attained at s=√(V⁴+V²Γ²).
This intentionally grants every pair the
maximum vacuum portal rate, even if its momenta
are far from the pole. It does not include hadronic
channels or non-vacuum propagator widths.

| Anchor and uniform cap | Scalar pair alone | Also grant the all-invariant three-lepton ceiling |
|---|---:|---:|
| 2 MeV; F≤1 | 13.16 | .2848 |
| 80 keV; F≤1 | 140.87 | 16.03 |
| 80 keV; F≤10 | 24.89 | 10.48 |
| 80 keV; F≤100 | .3460 | .3395 |

These are conditional lower bounds on the
stable contribution to Ωh² for the selected
reactions given the corresponding modeled start.
A uniform cap of at least about 170 would be
necessary for the 80 keV bound with leptons to
drop to .12; it would not be sufficient to reach
that value. The calculation neither derives
the actual occupation supremum nor constrains
new loss channels, nonstandard expansion, or
physical phase histories that depart from the
modeled anchor.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_post_2mev_elastic_reference.py EOS.late.dat feedback.json > late.json
python -W error warm_prompt_reheating_late_shape_free_floor_reference.py EOS.late.dat late.json
~~~

[Late histories and kick clocks](results/warm-prompt-reheating-post-2mev-elastic-2026-11-26.json.gz);
[shape-independent conditional floors](results/warm-prompt-reheating-late-shape-free-floor-2026-11-26.json).
