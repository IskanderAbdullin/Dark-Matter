# Refining the cold scalar shells changes sub-bin shape and local forward coefficients

The [two matched scalar-only histories](INDEPENDENT-EARLY-POLE-FINE-COLD-SCALAR-RESOLUTION-2026-11-28.md) have essentially identical total count at x=.05 on the same two-subcell charge grids. Aggregating their normalized distributions to the original 67 scalar shells gives 1.32% total variation. That comparison discards redistribution inside any original shell.

Treat the occupation as constant within each actual S2 or S4 shell. Then the number CDF is linear in q³ between successive shell edges. On the union of their edges, the CDF difference is also piecewise linear in q³. Its maximum occurs at one of these edges, and half the sum of absolute CDF-difference increments gives total variation of the reconstructed normalized counts. This yields **14.22%** TV on the union of all 145 shells and a maximum cumulative difference of **4.44 percentage points** at q≈3.625e−6. These are exact evaluations for the two declared histogram reconstructions, not uncertainties on a continuum solution. The [figure](results/early-pole-fine-cold-scalar-shape-2026-11-28.svg) and [receipt](results/early-pole-fine-cold-scalar-shape-union-2026-11-28.json) show the cold neighborhood and the all-shell metrics.

The [matched unrefined 4096-step control](INDEPENDENT-EARLY-POLE-FINE-NATIVE-SCALAR-STARTUP-CLOCK-2026-11-28.md) gives 34.70% union-shell TV from S2 and 33.24% from S4. The successive matched-grid comparisons (34.70%, then 14.22%) contract, but this is not an asymptotic convergence bound. Both are much larger than their respective common-native-bin distances; using a constant number density per unit q³ within each shell has exposed the redistribution hidden by coarse bins. Add the matched native state with `--matched-native` to reproduce these optional pairs in the receipt.

## Fixed transverse momentum, selected scalar momenta

Use each frozen state to reconstruct the charged distributions by their cell moments and solve the transverse pole at k=.0015. At five selected q, compare the forward coefficient with and without the scalar Bose multiplier. The ratio of the conditional coefficient is exactly the product of the raw kernel ratio and the multiplier ratio `(1+f_S4)/(1+f_S2)`.

| q | Conditional forward S4/S2 | Raw forward S4/S2 |
| --- | ---: | ---: |
| 2.6e−6 | 1.7712 | .9999158 |
| 3.5e−6 | 1.1748 | .9999170 |
| 4.0e−6 | 1.6327 | .9999173 |
| 4.7e−6 | .5851 | .9999181 |
| 8.2e−6 | .4610 | .9999240 |

At these points the change is chiefly the different scalar occupation, not changes to the raw transverse kernel. The conditional photon occupation that cancels the net source changes by only about 8–10 parts per million, since the same large scalar occupation multiplies both opposite process coefficients. This photon balance is a hypothetical value for an isolated mode pair, not a solved photon density. The [point receipt](results/early-pole-endpoint-scalar-refined-pole-points-2026-11-28.json) preserves each occupation and both quadrature comparisons (relative differences <6.5e−13). No q integral, physical collision rate, or convergence of either evolved history is inferred.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_endpoint_scalar_refined_pole_points.py \
  results/early-pole-fine-charge2-scalar{2,4}-x05-state-2026-11-28.json.gz
python -W error early_pole_endpoint_scalar_shape_figure.py \
  results/early-pole-fine-charge2-scalar{2,4}-x05-state-2026-11-28.json.gz /tmp/scalar-shape.svg \
  --receipt /tmp/scalar-shape.json
```
