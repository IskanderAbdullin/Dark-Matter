# Cold momentum transport and nodal bias measure different features of the same spectra

## Distance moved, not merely count reassigned

Compare the [time-matched endpoint histograms](INDEPENDENT-EARLY-POLE-FINE-NATIVE-SCALAR-STARTUP-CLOCK-2026-11-28.md) using the first Wasserstein distance W₁. For one-dimensional normalized counts this equals the integral of the absolute CDF difference in q ([Vallander, 1974](https://doi.org/10.1137/1118101)). Unlike total variation, it accounts for how far the count has to move. On the common union of shell edges each CDF is affine in q³. Integrate that exact form in q, splitting at its one possible zero inside each cell, and use reverse cumulative masses when both CDFs are close to one. Independent Gaussian integration and the signed mean identity check the result.

| Transition | W₁ for counts renormalized below q=1e−5 | W₁ / old cold mean | Cold TV | Absolute cold mean shift / W₁ |
| --- | ---: | ---: | ---: | ---: |
| Native to S2 | 2.59143e−7 | 6.715% | 34.70% | 81.36% |
| S2 to S4 | 3.87290e−8 | 1.062% | 14.22% | 68.82% |

The mean cold q moves down from 3.64840e−6 to 3.62175e−6 in the last step. The measured transport distance shrinks faster than TV over these two pairs. It is not an extrapolated spatial error bound. For a *common* smooth weight g on the normalized cold support, integration by parts yields `|ΔE[g]| ≤ max|g'| W₁`. This is not a bound for changing kernels, sharp reaction apertures, varying total number or a full transfer function.

## Higher moments depend on the representation, and velocity weighting is not q² weighting

Use the same fixed mass as the [spectral moment budget](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-SPECTRAL-MOMENT-BUDGET-2026-11-28.md). Treating each original cell weight as a point count at its geometric node reproduces total number. Compared with the exact piecewise-constant q²dq reconstruction, it undershoots the q² moment by 8.347% and the pressure numerator by 3.772–3.784% across all three cards. This is a representation difference at fixed saved occupations, not a quantified error in their dynamical evolution. Refining only the cold bins does not refine the warm grid that dominates those moments.

The number-weighted v² proxy, `q²/(q²+m²)`, saturates at high q. Consequently q≥.01 provides about 31.4–34.4% of this proxy, not the >99.8% obtained for pressure or the >99.9999% for q². In the normalized cold region, S2→S4 changes the proxy by −2.979e−7. Holding its mass fixed, the derivative bound on the cold interval gives an absolute change at most 1.016e−6 from the measured W₁. This verifies a specific smooth-weight example; it gives no CLASS or late-time velocity prediction.

The [transport receipt](results/early-pole-endpoint-scalar-shape-transport-2026-11-28.json) and [nodal/velocity receipt](results/early-pole-endpoint-scalar-nodal-moment-audit-2026-11-28.json) use the archived native 4096-step-start, S2, and S4 x=.05 states. From `kinetics/`, run the two corresponding scripts with those three gzipped state paths in that order.
