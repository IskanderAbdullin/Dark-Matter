# Same-flavor redistribution with separate energy conservation

The [screened unlike channel](INDEPENDENT-EARLY-POLE-SCREENED-RESPONSE-FEEDBACK-2026-11-26.md)
allows daughter/parent heat exchange. In contrast, each
same-flavor channel described here alters occupations without
changing either flavor's number or energy. They therefore
test a distinct way to change the narrow scalar production
window.

## Strict hard and soft-pole controls on disjoint domains

The hard subset uses the [massless Abelian Møller and Bhabha
amplitudes](https://indico.cern.ch/event/388801/contributions/1822106/attachments/1128409/1612259/HASCOLecturesDay1.pdf)
with their exchange and interference terms. Each pair initially
has `s≥1` and both vacuum fractions `|t|/s, |u|/s ≥ .15`.
The paired charge occupations are equal as in the archived
Dirac population. Relative to sixteen initial Dirac states:
like charge pairs have eight spin/charge states, with halves
for the identical incident pair and identical final pair;
opposite charge pairs have four states, with distinct charge
labels in both directions. There is also the common half
for the forward-minus-reverse integration over reciprocal
phase space. The finite dressed remap closes both energy
and an outward momentum triangle.

The second, separate, diagnostic uses the screened spacelike
t-exchange pole only at plasma-frame transfer
`1e-4 ≤ Q < .08`. For the like-charge process count both
permutation-equivalent t/u pole domains once, using the t
representative with twice its hard-sector population prefactor.
Add the Bhabha t piece with its distinct-charge prefactor.
The two contributions become .50 of the distinguishable
Dirac kernel on this t-oriented phase space. The soft and
hard tests are disjoint at the `s≥1` cut. They leave a gap,
and the soft pole calculation does not give the full
same-flavor interference or timelike terms. The response
in the pole sector is the frozen archived longitudinal and
transverse polarization.

At equal temperature within a given flavor, either test is
stationary independent of the other flavor's temperature
and chemical potential. The equilibrium fixture deliberately
sets different temperatures for parent and daughter.
Both channels also give an exactly conservative discrete
flavor energy and nonnegative entropy production up to
numerical accumulation.

## Instantaneous aperture response

The hard same-flavor contribution gives a positive summed
initial source tangent of about +45.2 (fine) or +3.9 (third)
at the coarse packet order. The much larger soft pole response
at the refined packet order is negative on both states:
−0.757 and −3.137, respectively, in units of `10^6`.
Most of this negative term comes from redistributing the
parent at unchanged parent energy. Photon and scalar states
are fixed for this directional derivative. Reducing the
soft transfer ceiling to .04 leaves its sign and magnitude
roughly unchanged in these finite comparisons.

At the same archived stage, adding this soft same-flavor
direction almost cancels the positive third-state tangent
from the unlike screened channel, but reinforces the negative
fine-state tangent. This is an instantaneous comparison of
linear directions, not a finite-time prediction.

## Shared implicit evolution

Group the three disjoint elastic packet graphs by final species.
The trial Jacobian is the sum of their Pauli derivatives;
the nonlinear residual recomputes each dressed map and volumes.
At twelve split steps and exposure `2e-6` on the native step
states, with photon pairs off, including the two same-flavor
soft pole sectors gives a relative scalar gain shift of
−3.17% (fine) and +11.25% (third) against the no-elastic control.
Unlike elastic alone gives about −9.6% and +9.9% at this
same step count. The initial tangent therefore does not
predict the relative effect of adding self scattering over
the whole finite interval. The calculations keep scalar
gain positive and close the energy and entropy ledgers.
No full temporal or collision-domain convergence is implied.

These are constrained local numerical experiments, not the
full elastic gauge rate or a late scalar result.

Run from `kinetics/` on one BLAS thread:

```sh
python early_pole_same_flavor_hard_redistribution_reference.py fine
python early_pole_same_flavor_screened_pole_reference.py third --level 1
python early_pole_same_flavor_screened_pole_reference.py third --cut .04
python early_pole_cold_screened_split_reference.py third --duration 2e-6 --count 12 --same-soft-strength 1
```

Raw: `results/early-pole-same-flavor-redistribution-2026-11-26.json.gz`.
