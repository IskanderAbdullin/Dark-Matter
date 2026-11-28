# Independent positive scalar volume nodes pass a short conservative live test

Extend the [frozen two-node quadrature](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-POSITIVE-CELL-QUADRATURE-2026-11-28.md) to an opt-in local collision chart with independent scalar occupations at the two positive volume nodes. There is no replacement of the normal driver. The charge nodes, geometry, and reaction logic are kept fixed. For an equal-frozen-input test, duplicate each selected cell's old constant occupation to its two nodes. Their normalized analytic weights are multiplied by the original stored volume; the largest adjustment between the stable factored exact volume and the legacy difference of cubes is 9.66e−14 relative. Cellwise count is preserved within floating-point precision.

From the fine S4 x=.05 state, compare the original 145 scalar nodes against using the two volume nodes in all 145 cells, or only in the 22 complete warm cells starting at q≥.01. The alternative full-node chart has 290 scalar nodes, while the warm-only one has 167. Existing conservation gates test eventwise energy deposition, paired rates and detailed balance, then advance 32 and 64 steps through the same short log-clock interval 5e−5.

| Changed scalar cells | 64-step scalar gain | Relative to single node | Maximum sampled 64-step continuity residual |
| --- | ---: | ---: | ---: |
| None | 6.31983129e−8 | 0 | 6.78e−17 |
| All | 6.31814047e−8 | −2.6754e−4 | 4.47e−17 |
| Warm only | 6.31982917e−8 | −3.3637e−7 | 1.45e−16 |

The local thermal detailed-balance relative defect is below 5.1e−15, and the imposed step-count change does not explain the reported gain contrasts. Pauli limits and number pairing pass at the tested states. The result is only a fixed inputs and short local evolution control, not evidence that a whole Gauss-node trajectory converges or that warm dynamics are physically negligible. Once the two nodal occupations are evolved independently, they are not a single piecewise-constant cell occupation.

The optional `early_pole_endpoint_scalar_gauss2_history.py` driver emits a distinct nodal checkpoint schema; its smoke test checks an empty startup, restart, and rejection of a resume with a changed node rule. It has not been substituted for other published trajectories. From `kinetics/`:
```sh
OPENBLAS_NUM_THREADS=1 python -W error early_pole_endpoint_scalar_gauss2_live_fixture.py \
  results/early-pole-fine-charge2-scalar4-x05-state-2026-11-28.json.gz
```
Add `--gauss-floor .01` for just the warm cells. The [all-cell](results/early-pole-endpoint-scalar-gauss2-live-all-2026-11-28.json) and [warm-only](results/early-pole-endpoint-scalar-gauss2-live-warm-2026-11-28.json) receipts preserve the charge and scalar counts and local diagnostics.
