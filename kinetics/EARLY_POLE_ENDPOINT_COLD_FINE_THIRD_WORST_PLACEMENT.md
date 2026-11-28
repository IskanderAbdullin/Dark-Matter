# How much does unmeasured within-panel location change the fine/third transport distance?

Keep the independently accepted normalized gross amount in every native panel. Conditional on those exact masses, the minimum interval-coupling cost is very different from the largest Wasserstein distance over all placements of those masses. Both are useful to separate the actual resolved ledger from unknown subpanel information.

For the upper value, use an ordered uniform probability coordinate. Each native panel covers a fixed subinterval of this coordinate. Its allowable quantile function is nondecreasing within the panel and bounded between the panel endpoints. With the other chart's quantile held fixed, write this function as a mixture of monotone two-level switches between endpoints. Convexity of absolute distance allows one such switch without decreasing the maximum. As the switch abscissa moves in the probability interval, the cost is convex, because the other chart's quantile is nondecreasing. Thus an extreme choice is constant at one of the two endpoints. Apply this to both charts. The probability overlaps form a forest because the ordered partition advances at least one new bin at every crossing. Binary dynamic programming maximizes the absolute displacement on this forest; exhaustive small tests verify the optimizer. Concentrating just inside the chosen endpoints approaches this upper value without requiring physical atoms.

| Gross response | Lower interval displacement | Sharp upper placement cost |
| --- | ---: | ---: |
| Forward number | 6.382744e−8 | 7.020727e−7 |
| Inverse number | 6.421399e−8 | 7.010814e−7 |
| Forward excess | 7.878841e−8 | 6.876201e−7 |

Even at fixed panel amounts, unresolved placement allows about eleven-fold variation in the forward distance. This is a limit of the retained kinetic data, not a prediction for the intrapanel profile or a physical momentum-transfer rate. The atom assignments are limiting mathematical witnesses and not a simulated distribution.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results; d=2026-11-28
python -W error early_pole_endpoint_cold_panel_transport_upper.py \
 --fine "$p/fine-d2-endpoint-continuum-cold"-{1,2}.jsonl.gz \
 --third "$p/third-d2-endpoint-continuum-cold"-{1,2}.jsonl.gz \
 --states "$p/early-pole-cell-endpoint"-{fine,third}-d2-x05-state-$d.json.gz \
 --archive "$p/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json" \
 --fixture "$p/early-pole-cold-fine-third-d2-native-cdf-$d.json"
```
