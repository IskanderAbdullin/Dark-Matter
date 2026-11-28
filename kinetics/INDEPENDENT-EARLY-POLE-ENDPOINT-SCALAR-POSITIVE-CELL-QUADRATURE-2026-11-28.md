# Two positive volume nodes recover frozen scalar moments without splitting the grid

The [single geometric-node audit](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-TRANSPORT-AND-NODAL-MOMENTS-2026-11-28.md) found an 8.35% second-moment and 3.77–3.78% pressure bias at otherwise fixed count. Here keep each actual saved chart, occupation and mass fixed and change only the quadrature on each shell. No collision step or scalar evolution is changed.

On a shell [a,b] write q=a+h u, h=b−a, with u in [0,1]. The shell volume measure is proportional to (a+h u)² du. Its normalized moments are obtained by integrating the three polynomial terms separately. Let μ₁, μ₂, μ₃ be its normalized moments of u and let s²=μ₂−μ₁². The symmetric two-by-two Jacobi matrix is
```text
[ μ₁       s                          ]
[ s        (μ₃ − 2 μ₁ μ₂ + μ₁³)/s²   ] .
```
The eigenvalues give two internal u nodes. The squares of the first components of the normalized eigenvectors are their positive normalized weights. Multiply them by the exact shell volume. This recovers the count and first three polynomial q moments of the assumed constant-in-shell occupation. The implementation verifies all four normalized u moments and ensures both nodes and weights remain interior and positive.

On the four matched native/S2/S4/S8 x=.05 histograms, the polynomial moments agree with independent 24- and 48-node reconstruction to floating-point rounding. The non-polynomial results at the same S2 mass are:

| Quantity | Relative bias, original geometric node (range across charts) | Largest absolute relative bias, two volume nodes |
| --- | ---: | ---: |
| Kinetic energy | −3.783% to −3.777% | 1.51e−10 |
| Pressure numerator | −3.784% to −3.771% | 1.22e−10 |
| Number-weighted v² integral | −5.707% to −.094% | 3.59e−8 |
| Inverse energy (Hartree proxy) | +1.38e−8 to +8.46e−7 | 3.39e−12 |

The inputs remain the same reconstructed histograms, including the fixed high-q part that dominates pressure. This is a cheaper way to evaluate smooth postprocessed moments on those inputs. It does not control within-shell occupation error, authorise replacing the collision chart, or prove a physical spectrum. An evolving method would separately need nonnegative occupation updates, exact paired budgets and a resolved comparison.

From `kinetics/`, run `early_pole_endpoint_scalar_positive_cell_quadrature.py` with the native 4096-step-start, S2, S4 and S8 gzipped x=.05 states in that order. The [receipt](results/early-pole-endpoint-scalar-positive-cell-quadrature-2026-11-28.json) reports all moment totals and observed differences against the higher-order references.
