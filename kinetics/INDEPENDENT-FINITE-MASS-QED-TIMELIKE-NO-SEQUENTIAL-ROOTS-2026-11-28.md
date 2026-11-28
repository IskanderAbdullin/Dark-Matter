# Timelike collective poles exclude bare charged-line sequential roots

The [spacelike photon test](INDEPENDENT-FINITE-MASS-QED-CUT-LOCAL-SEQUENTIAL-SUBTRACTION-2026-11-27.md) found occupied charged intermediate roots. A separate invariant argument rules them out globally on the [frozen timelike collective branch](INDEPENDENT-FINITE-MASS-QED-POLE-PAIR-APERTURE-2026-11-28.md). This is a kinematic separation, not a width or a rate.

Orient K to be future-directed with K²=μ²>0 and let P be a future on-shell external charged momentum with P²=m². In the rest frame of K, P·K=μE_P≥mμ. For the same bare mass on the internal leg the two signed denominators therefore satisfy

```text
(P+K)² − m² ≥ μ(2m+μ) > 0
(P−K)² − m² ≤ −μ(2m−μ) < 0,  provided μ < 2m.
```

The assertion covers all external magnitudes, both energy routings, and any angle. Equality in the dot-product bound occurs for comoving momenta. It is independent of occupations; it does not say that an on-shell external scalar reaction is otherwise supported.

In both frozen transverse runs, μ lies between the plasma gap and the maximum given by the transverse light-cone invariant. The upper bound is below even the light charged mass, so both magnitude polynomials are increasing throughout the representable interval. Evaluating at the lower plasma gap gives these uniform floors across both archived runs, in squared archive energy units:

| Charged leg | Lower bound for plus | Lower bound for magnitude of minus |
| --- | ---: | ---: |
| parent, m=.05 | 1.69154e−4 | 1.63617e−4 |
| daughter, m=.003 | 1.27515e−5 | 7.21473e−6 |

For the longitudinal timelike branch the pointwise inequalities hold at every μ>0. They give no positive uniform margin at its light-cone endpoint, where μ approaches zero. For example the comoving daughter minus denominator approaches zero along the checked β ladder. None of this argument holds in the spacelike sector, where the occupied sequential roots and their coherent coalescence remain separately documented.

The cloud replay independently checks a rationalized dot product constructed to avoid high-momentum collinear cancellation, its agreement with direct four-vector products on 3,456 sampled configurations, and saturation at comoving momenta. From `kinetics/`, run `python -W error finite_mass_qed_timelike_no_sequential_roots_reference.py` for the [JSON receipt](results/finite-mass-qed-timelike-no-sequential-roots-2026-11-28.json). These bounds assume the same external and intermediate **bare** mass; they do not choose a physical contour, interacting dispersion, pole statistical weight, or resolved-intermediate subtraction.
