# The aperture-local preparation attenuates the initial elastic source slope

Use the same held archived graphs, native hard mesh and initial cold
seed; replace only the two cold-active daughter cells. The initial
cold rate is 0.06312562923 rather than 0.07863619025. Differentiate the
cold source along the initial elastic moment direction, recalculating
dressed energies at both signs of three probe amplitudes. The raw
prepared snapshot independently verifies the exact initial local
hard and cold data.

| Initial held elastic derivative | Native | Local daughter |
|---|---:|---:|
| Two elastic channels (probe 3e-15) | -3.55691e6 | -4.59567e5 |
| Local minus native, unlike component | | +2.95383e6 |
| Local minus native, same component | | +1.43519e5 |

The local preparation lowers the initial rate but makes its elastic
direction much less negative. The difference in the two-channel
direction is +3.097e6 over all three probes (relative spread below
0.02%). About 95.4% of this initial differential is carried by the
unlike component. Restricting the derivative to the four originally
cold-active charged hard cells reproduces the differential within
0.01% over these probes. The outside daughter control differs from
native by only about -9e3 in the same derivative rather than the
large positive local shift.

These are directional derivatives in exposure units at the initial
held state, not slopes fitted to the full trajectory. Nonelastic
channels are not in this derivative. The result explains how a lower
initial rate can have a more positive elastic contribution to its
slope; it does not predict the time at which instantaneous or
integrated curves cross, or partition physical production.

From `kinetics/` use
`early_pole_daughter_initial_aperture_flux_reference.py` with the
paired n216 archive; `--local-result` additionally verifies a
matched local prepared split output. All three probes and the
mode-additivity diagnostics are in the paired JSON receipt.
