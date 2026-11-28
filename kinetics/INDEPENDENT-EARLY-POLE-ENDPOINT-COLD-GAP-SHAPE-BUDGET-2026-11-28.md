# The scalar Hartree gap is a weak probe of the reconstructed cold shape

The [matched S2 and S4 histories](INDEPENDENT-EARLY-POLE-FINE-COLD-SCALAR-RESOLUTION-2026-11-28.md) differ by 14.22% in reconstructed scalar histogram total variation at x=.05, yet yield almost the same scalar pole mass. Recompute the discrete Hartree source on both charts at the S2 pole mass m. Each scalar node of count N_i contributes N_i/E_i, where E_i=sqrt(m²+q_i²). Define cold shells by their upper edge q≤1e−5. In those shells,
```text
sum_cold N_i/E_i = N_cold/m − S,
S = sum_cold N_i q_i²/[m E_i(E_i+m)] .
```
The form with S avoids subtracting two nearly equal per-node inverse energies. Let bars denote S divided by the cold count. Then the charged-fixed change in cold source from S2 to S4 separates exactly as ΔN_cold·(1/m−bar_S4) + N_cold,S2·(bar_S2−bar_S4). The second term is the cold shape difference at the old cold count. Keep all shells above the cut as an independently summed contribution; do not expand them in q/m.

| Term in the S4 minus S2 gap source at the same m | Relative to the entire S2 source |
| --- | ---: |
| Cold count change | −9.3955e−7 |
| Cold shape, fixed old count | −2.6353e−10 |
| Above the cold cut | +4.9282e−9 |
| Total | −9.3488e−7 |

Cold nodes provide 99.99986% of the baseline inverse-energy source. The nonnegative deficit per particle in S is bounded between zero and its value at the cold cut. Even an arbitrary different distribution at fixed cold count within this support contributes at most 6.56e−5 of this baseline gap source through the shape term alone. The actual S2/S4 shape term is far smaller.

The gap residual is strictly increasing in m>0, with derivative at least 2m for this positive λ and nonnegative population. The measured source difference therefore bounds the relative shift of these two nodal roots by 4.68e−7; the observed shift is −3.1163e−7. This explains why agreement in the gap does not establish agreement in the cold spectrum. It does not control the separate [fixed-pole point response](INDEPENDENT-EARLY-POLE-FINE-COLD-SCALAR-RECONSTRUCTED-SHAPE-2026-11-28.md), for which local scalar Bose factors can differ greatly.

From `kinetics/`, run `early_pole_endpoint_cold_gap_shape_budget.py` on the two saved S2/S4 x=.05 states. The [receipt](results/early-pole-endpoint-cold-gap-shape-budget-2026-11-28.json) checks the gap equation and both accounting identities. The algebraic bounds are evaluated in floating point for this numerical gap; they are not physical rate or convergence bounds.
