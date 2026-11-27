# Add the momentum-triangle constraint to a manufactured scalar secant map

The [positive energy-only fixture](INDEPENDENT-EARLY-POLE-SCALAR-SECANT-POSITIVE-PROJECTION-2026-11-27.md) does not test whether its outgoing interpolation can represent the incoming total vector momentum. An optional diagnostic now retains the vacuum stencil's original incoming lab momentum magnitude P without changing its default packet data. In this manufactured variant every discrete outgoing candidate pair p₁,p₂ must satisfy |p₁−p₂|≤P≤p₁+p₂. These are the conditions for two vectors of those lengths to sum to the same incoming vector. The mixture of allowed pairs is adjusted to preserve the exact endpoint Hartree energy; if the allowed pair energies cannot bracket the target, that row is rejected.

| Scalar support | Retained Born weight, adjacent cells only | Retained Born weight, one added neighbor |
| --- | ---: | ---: |
| 1e−6 to 1e−5 | 97.88% | 99.49% |
| 1e−6 to 3e−3 | 93.91% | 99.09% |

These percentages refer to the unoriented positive weights of the initially accepted vacuum graph, not to a physical rate. The wider support's accepted row count is materially lower (5,534 of 8,392 on the adjacent-only construction). Including one neighbor substantially recovers weight but increases the broad fixture's entropy increment by about 26%. Recovering graph weights is not the same as converging the evolved distribution.

Every positively weighted outgoing pair meets the triangle condition, scalar count is unchanged at displayed precision, exact finite Hartree energy closes below 3.4e−21, and the independent entropy ledger agrees within 3.3e−12 relative. This constructs geometrically possible outgoing assignments; it does not recompute gauge or scalar matrix factors at the changed scattering angles, and energy is conserved in the weak interpolation mean rather than for each discrete candidate pair. The mobility is held in this fixture. The calculation is not yet a physically angle-matched collision operator and is not used in the live scalar-shape runs.

From `kinetics/`:

```sh
OPENBLAS_NUM_THREADS=1 python early_pole_scalar_triangle_secant_fixture_reference.py
```

Raw: `results/early-pole-scalar-triangle-secant-2026-11-27.json`.
