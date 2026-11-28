# Closed-reservoir depletion bound and the inherited stimulating band

These are two separate checks of the warm late pair channel.
Neither computes the fully coupled evolution.

## Energy and support bound with internal redistribution

The [selective-removal moment bound](INDEPENDENT-WARM-PAIR-SELECTIVE-LOSS-BOUND-2026-11-26.md)
assumed the survivors freely redshift. Relax that part here:
allow energy redistribution and inverse gains in an **isolated** active
dark reservoir with nonnegative pressure. Keep the incident
momentum cap (p\le P/b), the assumed photon floor, a future
scalar occupation cap F, and no external injection or re-entry
of the old tail. Only the stated pair channel may remove stable
number. New components produced within the active reservoir are
included in its energy ledger.

The total physical dark energy times (b^3) cannot increase under
these conditions. If x is the retained-and-new stable particle
count divided by its active initial count, the stable kinetic
energy comoving density is at most (E_0-mN_0x).
The incident second moment is at most this quantity times
(m+\sqrt{m^2+P^2/b^2}), by the exact identity relating p² and
kinetic energy. Use that inequality and the independent support
bound on the positive forward loss. Dropping the inverse contribution to the **loss
inequality** is safe here: it adds stable number,
while any stored kinetic energy from inverse gains is already
allowed by the weaker energy budget.

With (y=1-b^{-3}), integrate the smaller allowed loss rate.
The energy limit is proportional to
(x (E_0/(mN_0)-x)b^2[m+\sqrt{m^2+P^2/b^2}]).
The support limit is proportional to (x^2).
The ratio of the limits rises along the comparison trajectory
as b rises and x falls. Integrate to their first crossing,
then use the exact inverse-count support solution. A separate
no-event ODE calculation checks the crossing implementation.
The initial (E_0/(mN_0)) from the 1.206 MeV resolved state is 1.23923.

At nominal photon floor and F=1, equal charge populations give
a conditional final stable Ωh² floor **29.05**, versus **27.12**
for the support-only comparison at the same inert tail. The
loss-only fixed-moment result **52.10** is stronger but requires
the extra no-redistribution assumption. In the energy bound,
a cap must permit at least F=42.96 before the numerical lower
envelope reaches 0.120. Arbitrary isotropic charge abundances
and a 2% lower photon floor weaken this requirement to F=42.07.
These are necessary conditions in this restricted model, not
an estimate of actual stimulation, sufficient conditions for
the target, or exclusions of other histories.

## Occupation inherited at forward-pair final momenta

For incoming lab momenta below (P/b), the scalars' lab energies
must lie inside
`[hypot(m,P/b) - P/b, hypot(m,P/b) + P/b]`.
This follows by bounding the CM boost of each daughter by the
total parent momentum and then applying the incident cap. Use
both exact scalar mass-to-momentum transformations and b to
express the outer interval in the archive's momentum coordinate.
On that containing interval evaluate the exact maximum
occupation of the *existing reconstructed scalar only*, and
independently integrate its positive particle count.

| Anchor shape | All resolved scalar peak | Peak in the containing band at b=1 | Resolved scalar number in band |
| --- | ---: | ---: | ---: |
| 1.206 MeV | 7.056×10⁻⁴ | 1.758×10⁻⁷ | 2.304% |
| 1.000 MeV | 7.107×10⁻⁴ | 2.961×10⁻⁸ | 0.841% |
| 0.414 MeV, base | 7.056×10⁻⁴ | no overlap | 0% |
| 0.414 MeV, padded | 7.056×10⁻⁴ | 1.133×10⁻⁸ | 0.761% |

For example at the 1.206 MeV anchor, a global resolved peak by
itself overstates the maximum **inherited** occupation at possible
final momenta by about 4,000 times. The containing interval is
not an exact occupied support. The off-grid component and any
new pair products are not described by this evaluation; the
initial in-band value is **not** a cap on all future scalar
occupations in the energy argument above.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_pair_energy_headroom_bound_reference.py \
  "$D"/warm-coupled-pair-b10-28-shape-2026-11-26.json "$D"/warm-coupled-pair-b10-count-audit-2026-11-26.json
python -W error warm_pair_scalar_accessible_band_reference.py \
  "$D"/warm-coupled-pair-b10-28-shape-2026-11-26.json "$D"/warm-coupled-1MeV-12-and-20-shape-2026-11-26.json \
  "$D"/warm-late-paironly-{base,pad}-shape-2026-11-26.json
~~~
