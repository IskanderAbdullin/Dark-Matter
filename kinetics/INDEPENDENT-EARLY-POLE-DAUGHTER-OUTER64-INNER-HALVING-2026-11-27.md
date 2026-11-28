# Halving only the elastic step at 64 outer intervals improves the comparison

Use the same 64 held-response outer intervals, two local prepared
cells, cold seed, hard supports and total duration. Change each
elastic stage from one solve to two, with individual time widths
matched at exactly half the original, without changing where
slow half steps occur. The finest profile passes the original
terminal gates and the 5e-15 individual elastic budget
(maximum defect 1.78e-15).

| Profile | Mean local minus native rate |
|---|---:|
| 64×1 | +3.84570917e-4 |
| 64×2 | +4.39126558e-4 |
| Held simultaneous chart | +5.17245309e-4 |

The change +5.45556e-5 removes about 41.1% of the mean-rate
disagreement with the chart. The maximum sampled cumulative
gap decreases from 4.4392e-11 to 2.6123e-11. Cancelling a
hypothetical leading linear inner-width term gives
+4.93682199e-4, still 2.3563e-5 below the chart. That arithmetic
is a conditional diagnostic, not demonstrated time order.
The finite outer interleaving error need not vanish when the
inner step is refined.

Use the common subcycle tool with `--macro-steps 64 --elastic-substeps 2`
and separate native and local cases, then audit both profiles
and their simultaneous reference with
`early_pole_daughter_inner_halving_compare.py`.
