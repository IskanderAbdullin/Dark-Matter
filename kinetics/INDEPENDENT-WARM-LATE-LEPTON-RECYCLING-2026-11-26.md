# Allow the remaining vacuum charged-lepton portal after 10.442 MeV

The earlier [energy envelope](INDEPENDENT-WARM-ISOLATED-RECYCLING-ENVELOPE-2026-11-26.md)
and [known-tail extension](INDEPENDENT-WARM-KNOWN-EXTERIOR-RECYCLING-2026-11-26.md)
switch off new visible sources at the archived anchor. This
test allows the specified positive **vacuum charged-lepton
pair injection and elastic heating** after it. The visible
entropy history remains fixed. Dark momentum redistribution,
tail reentry, inverse decays and repeated recycling can occur
in the same restricted internal network. Ignore energy and
number losses to the visible sector to get an upper bound.

## Bound the supplied energy rather than select a new spectrum

Use the [project's energy-source identity and heat ceiling](https://github.com/IskanderAbdullin/Dark-Matter/blob/90c141722376dab64e4cd10422d545474817787c/kinetics/PHI-BOSE-PORTAL-RECYCLING-ENERGY-BOUND-2026-11-24.md).
For forward pair injection, Maxwell–Boltzmann initial weights
and no final blocking upper-bound the stated FD rate. After
writing u=√s/T, the energy moment uses
`u⁶[2K₁(u)+uK₀(u)]`, with exact full integral 3072.
The massive flux and phase factors are each bounded by 1,
so the six lepton/species combinations share one massless
ceiling. Below √s=1 GeV bound the vacuum propagator by
`(m_V²−1 GeV²)⁻²`. Above that split allow the entire
Breit–Wigner peak through its width floor, but bound the
thermal tail via `K_ν(u) ≤ exp(−u) sqrt[π/(2(u−1))]`
for ν=0,1. This follows directly from the cosh representation
and leaves only positive upper incomplete gamma terms. The
tail correction is below 3.1×10⁻²⁶ of the contact part.
As a normalization diagnostic the direct massive MB fold at
the anchor is 13.74% of this deliberately loose bound.

For elastic scattering the same vacuum model gives
`Q_heat ≤ 288 C T⁵ ρ_D/(π³ m_V⁴)`, with `C` the squared
vector coupling product. This uses the second energy moment
of the six charged species with two spin states each. In
every visible-EOS interval, extrema of entropy and radiation
coefficients give analytic upper integrals of the required
temperature powers over expansion. After b=500 use only the
photon expansion floor and at least two visible entropy
degrees of freedom. Integrate the pair inflow and heating over
all future time; credit all resulting energy as available
already at the anchor when bounding internal recycling.
Each injected heavy can also decay once; direct net particle
input including these daughters is at most pair energy divided
by the stable mass. This keeps direct injection distinct
from later heavy formation in the internal ledger.

| Anchor | Positive pair energy / initial core energy, upper | Elastic integrating-factor exponent, upper | Scalar Ωh² upper if total is .120 |
| --- | ---: | ---: | ---: |
| q16 | .020514 | 7.266×10⁻⁵ | .00292489 |
| q32 | .020502 | 7.266×10⁻⁵ | .00292398 |

Even allowing these vacuum lepton inputs, the scalar fraction
is at most 2.44% if the total is .120. This does not predict
the total, include late hadrons, include new scalar
number-changing processes, or cover energy transfer before
the anchor that was absent from the archive. The analytic
massless ceiling and pole-tail bounds are not a replacement
for an actual coupled source calculation.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_late_visible_recycling_envelope_reference.py EOS.dat \
  "$D"/warm-isolated-recycling-energy-envelope-2026-11-26.json \
  "$D"/warm-known-exterior-recycling-envelope-2026-11-26.json
~~~
