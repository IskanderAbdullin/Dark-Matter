# Frozen hard Compton heat on imposed Bose photons

The [screened interspecies elastic study](INDEPENDENT-EARLY-POLE-SCREENED-RESPONSE-FEEDBACK-2026-11-26.md)
has no Compton exchange. Here independently add a *frozen*
finite vacuum-hard diagnostic of `γ + f ↔ γ + f`.
Use the original two archived fermion distributions and an
imposed Bose photon occupation. Include full stimulated emission
and Pauli blocking in both directions. The distinct initial
states have two photon and four Dirac states; the extra
factor one half compensates for applying forward-minus-reverse
over the complete reciprocal domain.

Select `s ≥ 1` and a central CM angular gate
`τ=(1+cosθ)/2 ≥ .15` where θ is the photon scattering angle.
In the massless case `τ=−u/s`.
The [spin-summed vacuum amplitudes](https://fma.if.usp.br/~burdman/QFT1/lecture_18.pdf)
give the massless weight `(1/τ+τ)/(16πs)` for the differential in
`cosθ`, with the common `g⁴` stripped. The finite-mass repeat uses
the same angle cut, not the same exact invariant cut. Its spin
normalization is checked against the rest-frame Klein–Nishina
angular bracket. The angular integral before final radial cuts
also agrees with its independent full-angular normalization.

Both incoming and outgoing momenta must lie in [0.45, 8].
The outgoing azimuth quadrature is split at every relevant
archive fermion edge and at the photon box boundaries.
A common-temperature test with independent chemical potentials
balances to roundoff. A hotter photon bath gives positive net
fermion heating on that same phase-space window.

At imposed photon temperature 1 (in reference units), the
medium quadrature gives signed heat into each fermion, in
`T_ref⁴` per reference log-scale interval at `g=.005`:

| Source spectra | Parent | Daughter |
| --- | ---: | ---: |
| Fine | −15.72 | +12.11 |
| Third | −16.86 | +11.87 |

The finite parent mass moves the same coarse value by less than
0.7%, and the smaller daughter mass shifts it by less than 0.1%.
Neither changes the sign. Lowering the angular floor to .10
alone changes the fine pair of values to −20.10 and +15.61;
widening only the radial box to [.30, 10] changes them to
−17.85 and +13.49. The chosen central window is not the whole
scattering phase space.

Changing the imposed photon temperature is more important
than the finite-mass correction in this window. At .65 the
daughter heat is still about +1.5, while the parent cools.
At 1.075 both archives have negative parent heat, but at
1.10 both have positive parent heat. Daughter heat stays
positive on the sampled .65–1.20 temperatures. These brackets
describe the imposed-photon comparison, not the actual photon
history.

These are **instantaneous exchanges against an imposed photon
state**. The photons are not being evolved and no self-consistent
photon temperature is measured. In a live model the opposite
energy change goes to photons, and neither species count can
change in this channel. The heating estimates cannot be turned
into a cold-scalar rate by addition to the unlike-fermion heat:
the final state and its narrow aperture response must be
computed with a common energy and occupation update. Soft,
collinear and other gauge contributions remain absent here.

Run from `kinetics/` (one BLAS thread):

```sh
python early_pole_hard_compton_heat_reference.py --level 1
python early_pole_hard_compton_heat_reference.py --flavor parent
python early_pole_hard_compton_heat_reference.py --temperature 1.05 --without-fixtures
python early_pole_hard_compton_heat_reference.py --bounds .3 10 --angular-floor .1 --without-fixtures
```

Raw: `results/early-pole-hard-compton-heat-2026-11-26.json.gz`.
