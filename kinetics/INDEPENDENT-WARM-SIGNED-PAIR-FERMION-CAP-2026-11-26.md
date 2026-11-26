# Test the warm fermion floor without imposing a scalar occupancy cap

The [previous floor](INDEPENDENT-WARM-FERMION-SURVIVAL-WITH-LEPTONS-2026-11-26.md)
requires a uniform bound on final scalar Bose occupancy.
For a spatially homogeneous isotropic gas, instead bound
**fermion occupations per spin**, for both masses and both
charges, by a common `ε≤1/2`. Allow a scalar distribution
with arbitrary occupancy and finite energy. This is a
different conditional comparison, not a proof that either
required occupancy condition is realized.

## The signed operator removes the quadratic Bose term

At matched phase space, the signed forward-minus-inverse
population factor is exactly

`f₁ f₂ (1+fₐ+f_b) + (f₁+f₂−1) fₐ f_b`,

where 1 and 2 are the fermions and a,b the scalars.
The last term is nonpositive under the fermion cap and
cannot increase net fermion loss. The first term has a
vacuum part, treated as before, plus two terms linear in
actual scalar occupancy. This argument uses **signed**
loss: omitting the inverse before expanding would leave
an unbounded positive quadratic term.

## Bound the linear terms by energy and fermion occupation

For initial fermion mass d and exchanged fermion mass e,
the same spin-trace ceiling gives `|M|²≤2 y⁴ s/e²≤8 y⁴ E₁ E₂/e²`.
For one outgoing scalar of energy E, integrate out its
partner with momentum conservation. The angle of the
second fermion gives a delta factor bounded by
`π/[p₂ |p₁−k|]`. Isotropic averaging of the first
fermion gives `4π/max(p₁,k) ≤4π/p₁`. With the two terms,
four initial spin states and identical outgoing scalar
factor, the resulting event bound is
`(y⁴/π³e²) ∫ dΠ_k f_k J₁ J₂`, before extra threshold
restrictions, with `J_i=∫ p f_i dp` and
`dΠ_k=d³k/[(2π)³ 2E]`.

Split at `E=d/(2η)` for η>1. Below it, the partner
scalar needs energy at least `d²/E`; at least one parent
therefore has energy L≥`d²/(2E)≥ηd`. Its radial
tail is bounded by `J_tail≤π²ρ_charge/[L² sqrt(1−η⁻²)]`.
The E² cancels the soft scalar denominator against
the scalar energy measure. Above the split, use
`dΠ_k≤2η²/d²` times that energy measure.
A filling argument under the per-spin cap gives
`J≤ε^(1/3)(3π² n_charge)^(2/3)/2`, attained by
putting the maximal allowed occupation at lowest p.

Combine charges and species conservatively. In the soft
term use `ρ_f ρ_φ≤ρ_dark²/4`; in the hard term use
`ρ_φ≤ρ_dark` and the equal-charge maximum of the radial
product at fixed number. The two moment integrals have
weights `b⁻⁵/H` and `b⁻⁴/H`. Work on the last descent
below the conservative initial core fermion count, on
which that count upper-bounds both radial integrals.
Reserve the integrated additional linear losses as
additive debits after the usual vacuum quadratic
comparison, just as for the high-invariant lepton loss.
Positive pair injection, elastic heating, and all
known tails use the earlier energy budget.

| q32 uniform fermion per-spin cap ε | Lower fermion Ωh² with no scalar cap |
| --- | ---: |
| .0001 | 21.05 |
| .001 | 13.09 |
| .003 | 3.83 |
| .01 | No positive bound from this comparison |

The q32 conditional exclusion of total .120 persists
up to a cap of about .00405. The archived retained
light-fermion peak at the 10.442 MeV q32 anchor is
about 5.88×10⁻⁵, but this does **not** bound later light
or heavy occupation or the omitted tail. Unresolved
momentum anisotropy, new channels and hadronic losses
remain outside the comparison. The zero-cap entry in
the numeric artifact is only a formal check against
the nonstimulated result, not a possible nonempty gas.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_signed_pair_fermion_cap_reference.py \
 "$D"/warm-isolated-fermion-survival-bound-2026-11-26.json \
 "$D"/warm-fermion-survival-with-leptons-2026-11-26.json \
 "$D"/warm-isolated-recycling-energy-envelope-2026-11-26.json
~~~
