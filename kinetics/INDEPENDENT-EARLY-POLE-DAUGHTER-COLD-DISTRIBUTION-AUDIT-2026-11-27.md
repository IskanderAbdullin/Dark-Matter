# The cold production increment and the populated terminal spectrum respond at different scales

Use the accepted native/local held-chart mode histories
separately on the eight- and 32-bin grids. Each pair has exactly
the same initial counts, node positions and integration settings.
Both individual production increments are nonnegative in all
reported bins, so both the increment and the terminal cold
distribution admit separate probability normalizations here.

| Cold bins | Local excess relative to native increment | Local excess relative to native terminal count | TV of normalized increments | TV of normalized terminal spectra |
|---|---:|---:|---:|---:|
| 8 | +0.7778% | +9.95e-9 | 2.230e-4 | 2.871e-10 |
| 32 | +0.7642% | +9.77e-9 | 1.168e-4 | 1.503e-10 |

TV is half the sum of absolute differences of per-mode probabilities
on the identical nodes in each pair. The numerical implementation
normalizes the paired difference in cancellation-safe form. On
the 32-bin grid the relative mean momentum shift is +41.54 ppm
for the increment, versus +5.33e-11 for the populated terminal
spectrum. About 80.55% of the additional count lies in the two
original groups numbered five and six. The receipt also measures
the one-dimensional normalized transport distance on each grid.

These are different denominators over the same short model interval.
An approximately 0.76% increase in the newly added cold count is
not a 0.76% change of the populated spectrum. These results do
not estimate a relic abundance, later transfer or a long-time
error bound. Use `early_pole_daughter_cold_distribution_audit.py`
with the accepted native/local mode histories at each grid.
