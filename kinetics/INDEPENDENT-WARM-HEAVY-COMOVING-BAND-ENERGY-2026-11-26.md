# Bound heavy energy in a finite comoving band during the warm evolution

The [direct light source test](INDEPENDENT-WARM-LIGHT-DIRECT-PAIR-DOSE-2026-11-26.md)
would benefit from a justified heavy energy allowance much
smaller than the whole dark-sector budget. Keep the assumed
light per-spin cap δ, the same finite gross dark-energy
allowance (E_A), and the [bound on other heavy pair inputs](INDEPENDENT-WARM-HEAVY-OCCUPATION-BARRIER-2026-11-26.md).
Assume the heavy band is initially empty and no unbounded
elastic transfer moves heavy fermions into it or supplies
them energy. Allow arbitrary scalar occupations with finite
energy. Bound only the heavy modes below a declared
comoving cut (P_c), normalized to the 10.442 MeV anchor.

A simple mode-by-mode bound uses the exact smallest scalar
daughter energy (e_-(P)) and the lab angle Jacobian
(M/(2 p P_*)). Along a characteristic (p=P/b) the
physical scalar energy-density ceiling falls as (b^{-3});
that factor more than cancels the (b/P) angle Jacobian.
Conservative shell maxima give a finite band energy bound,
but it is too loose if summed over the whole interesting range.

The aggregate comparison instead integrates the positive
inverse source before bounding it. Its energy source divided
by (Gamma M), integrated over both charges and spins, is

`M/(π² P_*) ∫ dE_l de (E_l+e) f_l(E_l) f_φ(e) 1_kin`.

The (E_l+e) factor comes from the energy weight of the
tagged heavy particle. The on-shell delta sets the relative
angle of the incoming pair. A separate coordinate integration
of the pointwise source on three smooth inputs agrees with
this aggregate normalization at better than (2×10^{-15})
relative. Retain the exact lower light energy from the
two-body constraint at fixed e and the upper requirement
(E_l+e≤sqrt(M²+P_c²)).

For fixed scalar energy, the ratio of positive light-source
weight (E_l+e) to the light energy measure (p_l E_l²) is
decreasing in (E_l). The worst allowed occupation therefore
fills the lowest accessible light energies first under δ.
Bound the filling endpoint from its energy allowance.
The arbitrary scalar occupancy can spend the remaining energy
in the worst allowed cell. Enclose both energies and the
density split by positive endpoints so the same full energy
cannot be assigned independently to both species. The input
upper at the anchor also bounds all later comoving source
weights. For the band energy Y, decay supplies loss at least
(Gamma M(1−δ)Y/sqrt(M²+P_c²)), while expansion lowers Y.
This bounds Y by the maximum input divided by this minimum
relaxation coefficient for an initially empty band.

| q32 with δ=1/2, comoving cut | Fraction of (E_A) allowed in heavy modes **below** the cut |
| --- | ---: |
| .10 GeV | ≤.00013587 |
| .30 GeV | ≤.0013711 |
| .60 GeV | ≤.0052538 |
| 1.00 GeV | ≤.0281222 |
| 1.40 GeV | ≤.0957330 |

At 1 GeV the below-cut bound is 2.8123% of the budget.
To meet the earlier partial direct-light diagnostic of
4.8336%, one would still need an independent bound of
about 2.021% of the budget **above** that cut. At .60 GeV
the corresponding above-cut allowance is about 4.308%.
The high band is not zeroed or claimed absent. These
conditions do not include light decay/elastic
redistribution and are not proof of a light occupation cap.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_heavy_low_comoving_energy_bound_reference.py \
 "$D"/warm-signed-fermion-with-charged-pions-2026-11-26.json \
 "$D"/warm-heavy-internal-source-barrier-2026-11-26.json \
 "$D"/warm-light-direct-pair-dose-2026-11-26.json > /tmp/heavy-pointwise-band.json
python -W error warm_heavy_band_aggregate_energy_bound_reference.py \
 "$D"/warm-signed-fermion-with-charged-pions-2026-11-26.json \
 "$D"/warm-heavy-internal-source-barrier-2026-11-26.json \
 /tmp/heavy-pointwise-band.json
~~~
