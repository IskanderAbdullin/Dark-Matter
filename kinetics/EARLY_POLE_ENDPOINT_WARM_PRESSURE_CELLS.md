# Native-cell attribution of the warm pressure change

Use the same four S4 and W2 endpoints as the [warm observable quantiles](EARLY_POLE_ENDPOINT_WARM_OBSERVABLE_QUANTILES.md). Integrate the evolved piecewise-constant spectra at the S4 mass, then collect subcell moments on the unchanged native scalar cells. Write each native contribution as count times mean pressure weight. For count `N` and mean `w`, the exact difference is `delta(Nw) = delta(N) * mean(w) + mean(N) * delta(w)`, where means average the two charts. This symmetric bookkeeping includes the count/shape interaction without assigning it all to one chart. An empty cell uses the populated chart's mean for both, since its own mean is undefined. Nothing here assigns a pressure contribution to a time integral or a reaction channel.

At x=.05, the total numerator change of -6.23477% divides into -3.04112 percentage points from changing counts between native cells and -3.19366 points from changing the pressure weight within those cells. Thus a match in the total number does not fix either term. The three adjacent q cells spanning [1.31278,3.96368] give -4.09296 points (65.65% of the signed change). The net deficit stays nearly the same over the four sampled endpoints, split between -3.04 to -3.13 points from cell counts and -3.19 points from within-cell shape. These contributions use reconstructed cell distributions, not collision-source attribution or a physical pressure forecast.

With the endpoint harness dependencies available, run from `kinetics/`:
```sh
python -W error -m unittest test_early_pole_endpoint_warm_pressure_cells
for i in 01 02 04 05; do
  python -W error early_pole_endpoint_warm_pressure_cells.py \
    results/early-pole-fine-charge2-scalar4-x${i}-state-2026-11-28.json.gz \
    results/early-pole-fine-charge2-scalar4-warm2-x${i}-state-2026-11-28.json.gz --x 0.${i}
done
```

The stronger integration compares 24 and 48 nodes relative to the gross contribution in each cell. Count must survive rebinning to the native edges, values must close under the symmetric difference identity, and the true unchanged cells may not show a change in within-cell shape. The four totals also independently agree with the pressure integral from the weighted quantile check.
