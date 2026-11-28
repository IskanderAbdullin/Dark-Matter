# Locate the warm collision-quadrature bias

At fine S4 x=.01 and .05, hold the scalar histogram, charged distribution, scalar mass and clock fixed. Reduce the warm pressure collision source into the 22 whole native cells above q=.01. Compare geometric W1, W2 and W4, positive two-node volume Gauss, and 32- and 64-point volume rules, with a paired 4/8 charge radial check. Forward, inverse and signed coefficients close within each bin, and the high-order sum closes against the independent whole-support reduction.

At x=.05 the three adjacent cells q=1.31278 to 3.96368 carry 57.8% of the 64-node net source and give 66.9% of W1's signed bias. The same cells account for 65.65% of the independently measured change between evolved S4 and W2 pressure snapshots. This overlap is consistent with a momentum-localized numerical effect; it doesn't establish a decomposition of the time evolution.

A frozen selective rule ranks cells solely by the absolute W1 discrepancy at x=.01, then uses the same mask for x=.05. Replacing the top eight W1 cell contributions with Gauss reduces the held-out signed error from +9.221% to +0.267%; the sum of absolute cell errors is still 0.472% of net reference. At 12 cells those errors are +0.123% and 0.250%; full Gauss gives +0.082% and 0.205%. This is a diagnostic that combines already computed source contributions. It neither changes the time-evolution chart nor establishes a recommended production mask. In particular, an early small signed error at 12 cells partly cancels errors and does not transfer unchanged to the later state.

The unsigned sum of 32-to-64 discrepancies is under 1.20e-5 of total net pressure source at each time (finite-mesh comparison, not an extrapolated bound). A uniform cell-accuracy gate of 4e-5 relative to gross source **fails** both times, with the largest defect about 1.04e-3 in q=[17.30,25], which contributes only about 6.3e-7 of the late reference net pressure. Thus the available reference is adequate for this aggregate comparison, not for a uniformly accurate fine-grained high-q source.

With the endpoint dependencies available from `kinetics/`:
```sh
python -W error -m unittest test_early_pole_endpoint_warm_source_cell_hybrid
for i in 01 05; do
  python -W error early_pole_endpoint_warm_source_cells.py \
    results/early-pole-fine-charge2-scalar4-x${i}-state-2026-11-28.json.gz \
    > results/early-pole-endpoint-warm-source-cells-x${i}-2026-11-28.json
done
python -W error early_pole_endpoint_warm_source_cell_hybrid.py \
  results/early-pole-endpoint-warm-source-cells-x01-2026-11-28.json \
  results/early-pole-endpoint-warm-source-cells-x05-2026-11-28.json
```

As an independent check of the frozen-mask construction, `early_pole_endpoint_warm_live_mask.py` builds the actual mixed 153-node chart (native W1: 145, full Gauss: 167). Its paired 4/8 radial collision evaluations at both snapshots reproduce the hybrid’s native-cell affine ledgers exactly in these runs, while conserving energy and retaining each native scalar cell count on lift. This verifies the mixed chart for the frozen test, without executing an adaptive or masked history. For reproduction, pass the corresponding S4 state, the cell-source JSON for that time and the hybrid JSON to the script, with its default `--cells 8`.
