# The local held-chart shift persists across eight and 32 cold bins

Hold the original hard supports, response, two elastic channels
and the exact two-cell local daughter substitution. Increase
only the cold grid from eight bins to 32 nested bins over the
same momentum aperture. The four-way regrouping of fine seed
counts differs from the original coarse seed by at most
8.67e-19. The reported simultaneous histories use the same
time points and tolerances, and all four states pass the
original ledgers.

| Cold bins | Native mean | Local mean | Local minus native |
|---|---:|---:|---:|
| 8 | 0.06649893898 | 0.06701618429 | +5.17245308e-4 |
| 32 | 0.06650019808 | 0.06700839244 | +5.08194357e-4 |

The local mean shift falls by 1.75% under this refinement but
remains positive. Every individual refined bin has a positive
terminal local-minus-native count, as does every original bin.
The largest regrouped difference is in original bin six, so
matching total results does not make the two spectra identical.
Figures in this table use recorded per-mode count differences;
summation order can change the last digits of total rates.

Use the matched 32-bin native/continuous split archives as
the two branches of one bundle. Run the simultaneous chart
with `--record-mode-history` and either native defaults or
`--local-from-paired-archive`; the latter keeps the exact two
original daughter cells. Compare with the eight-bin
native/local histories using
`early_pole_unsplit_local_cold_grid_compare.py`.
These are numerical controls within one held model, not a
continuum error bound or physical abundance claim.
