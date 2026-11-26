# Can late recycling overcome the count deficit with redistribution?

This is an archived-energy continuation of the project's
[recycling energy bound](https://github.com/IskanderAbdullin/Dark-Matter/blob/90c141722376dab64e4cd10422d545474817787c/kinetics/PHI-BOSE-PORTAL-RECYCLING-ENERGY-BOUND-2026-11-24.md).
It tests a different question from the
[passive redshift integral](INDEPENDENT-WARM-BALLISTIC-HEAVY-EXPOSURE-2026-11-26.md):
start at the matched warm 10.442 MeV snapshots, allow arbitrary
dark momentum redistribution and stimulation within the specified
reaction network, and cease all external injection, heating and
reentry. Use only the represented energy as the initial isolated
reservoir. The particles in this reservoir are stable χ₂ and φ.
Positive pressure implies `ρ_D(b) ≤ ρ_D(1)/b³` for the combined
reservoir, even if pair reactions or prompt/inverse decays occur.

## Energy-weighted rate and count ledger

1. For χ₁ anti-χ₁ → φφ, a single t-channel spin trace is
   `(m₂²−t)s/2 − (m₁²+t−mφ²+2m₁m₂)²/2`; likewise for u.
   The transfers t and u are nonpositive in the center of mass.
   Cauchy bounds the squared t/u sum by twice the sum of squares,
   so the angle-integrated vacuum `σ v_CM` is at most
   `βφ y⁴/(8π m₂²)`, keeping the outgoing velocity factor.
   For a tighter integral, write each virtual denominator as
   `A ± D cos θ`, with
   `A=m₂²+s/2−m₁²−mφ²` and
   `D=2 sqrt[(s/4−m₁²)(s/4−mφ²)]`.
   Its reciprocal integrates to `2 atanh(D/A)/D`, at most
   `2/sqrt(A²−D²)` by the derivative comparison
   `atanh z ≤ z/sqrt(1−z²)` for 0≤z<1.
   The squared geometric denominator is
   `m₂⁴ + m₂²(s−2m₁²−2mφ²) + (m₁²−mφ²)²`, increasing
   in s and minimized by `(m₂²+m₁²−mφ²)²` at threshold.
   Replace `m₂²` by `m₂²+m₁²−mφ²` in the forward bound.
2. Detailed balance contributes `8 β₁/βφ` for the inverse scalar
   channel, before blocking. The energy-weighted lab rate obeys
   `σ v_lab/(E E′) ≤ (4/s) σ v_CM`, since
   `v_lab/v_CM = s/(4 E E′)` and `E E′ ≥ s/4` for the identical
   incoming masses. The maximum of `(4 m₁²/s) β₁` is
   `2/(3√3)` at `s=6 m₁²`. Thus the scalar number-increase rate
   from this route is at most `Kφ ρφ²/m₁²`, with
   `Kφ=(2/3√3)y⁴/[π(m₂²+m₁²−mφ²)]`: the 1/2 for identical initial scalars
   cancels the two heavy parents per event.
3. For endothermic light-fermion pair conversion, use the
   upstream t-channel ceiling, identical initial-channel weights
   and the same-charge Cauchy bound. The combined stable-channel
   coefficient for number increase is at most `2 KD` in front
   of `ρχ₂²/m₁²`. The two energies are disjoint, so the total
   coefficient multiplying `ρ_D²/m₁²` is
   `max(Kφ, 2 KD) = 1.31980×10⁻¹¹ GeV⁻²`.
4. Every net decay beyond an inverse decay consumes a heavy
   parent. With none initially, integrated net decays cannot
   exceed the integrated number of gross newly formed heavy
   parents. Any heavier-pair loss only reduces this number.
   Consequently this bound does not assume one pass, small
   scalar occupations, or preserved momentum support.

Integrate `b⁻³/H` in `d ln b`. The released central
[EOS table](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat)
provides the visible clock through b=500. Bound every interpolated
temperature interval by its lower radiation-energy and higher
entropy endpoint, using the no-overshoot PCHIP interpolation.
Use the weaker photon-only `H ≥ Hγ(1)b⁻²` floor afterward.
Additional nonnegative energy increases H at fixed visible
entropy history. Independent temperature-coordinate quadrature
checks, but does not replace, the interval envelope. A separate
gamma-matrix calculation checks the single-channel spin trace
to 2.35×10⁻¹⁵ in units of `(m₂²−t)s` over 510 probes.

| Input at 10.442 MeV | Upper final count / initial count | Upper Ωφh² if total is .120 | Count ratio needed for half |
| --- | ---: | ---: | ---: |
| q16 | 5.347 | .002807 | 114.29 |
| q32 | 5.352 | .002809 | 114.29 |
| q32, +3% initial energy and −2% H throughout | 5.711 | .002998 | 114.29 |
| q32, +5% initial energy and −2% H throughout | 5.896 | .003095 | 114.29 |

The count-to-Ω conversion uses the somewhat larger initial
inclusive count normalization from the
[closed-count audit](INDEPENDENT-PORTAL-CLOSED-NUMBER-BUDGET-2026-11-26.md);
it does not grant the missing particles or their energies reentry.
At the nominal input, even this loose upper envelope leaves φ
below half of a total .120 *if that total is reached*. It does not
predict that the total reaches .120. The nominal q32 scalar
upper is 2.35% of that total. The original uniform-angle bound
was much looser (about 102 initial counts) and its joint +5%
energy/−2% H variant crossed the 114.29 threshold; the analytic
angular bound removes that sensitivity in these controls.
Neither the bound nor the passive
calculation constrains unknown earlier heating, later inflow,
other number-changing vertices, or the full physical model.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_isolated_recycling_energy_envelope_reference.py EOS.dat \
 "$D"/warm-and-low-closed-number-budget-2026-11-26.json \
 "$D"/warm-q{16,32}-10p442-heavy-input-shape-2026-11-26.json
~~~
