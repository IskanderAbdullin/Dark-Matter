# Evolve the original hard aperture with shared angular-aligned cold modes

The 78-node [initialization](INDEPENDENT-EARLY-POLE-DAUGHTER-ANGULAR-MODE-INITIALIZER-2026-11-27.md)
is now propagated, with both paired hard-state interventions,
on the original hard depth-one geometry and the nominal
unsplit first-window tolerances. The cold nodes, weights,
initial numbers and initial angular panel hash agree between
arms. Cold counts remain positive and reproduce each reported
scalar count exactly in the original summation order. An
independent compensated sum differs by at most six ulps.

The final paired cold mean is +5.0765147294e-4, versus
+5.0819435118e-4 on 32 logarithmic volume-centroid modes
(a reduction of about 0.107%). All 78 individual terminal
paired increments are positive. The cumulative sign change
remains bracketed by [1.250, 1.458]e-7, as on the ordinary
mesh. Maximum energy closure error in the angular pair
is 8.89e-16. The earlier shared static quadrature error is
not an error bar for this trajectory. In particular this
result does not yet compare angular evolution on the more
refined hard aperture.

Use `early_pole_daughter_angular_held_chart.py` with the
original paired initializer and shared angular mesh receipt
for both arms, adding `--local` for the daughter intervention.
Pass the two angular outputs, the previous two nominal
cold-32 outputs and the 78-node initializer receipt to
`early_pole_daughter_angular_held_chart_compare.py`.
Both accepted compressed angular raw histories are retained.
