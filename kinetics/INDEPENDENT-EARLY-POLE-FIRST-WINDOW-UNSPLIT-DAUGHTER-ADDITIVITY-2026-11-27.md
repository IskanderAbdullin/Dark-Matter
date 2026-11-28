# The local preparation offsets most of the outside daughter effect in the unsplit test

The two hybrid inputs changed disjoint sets of daughter moments on
common native support. Exactly combining their initial arrays
reproduces the entire continuous daughter substitution with the
native parent and photon preparation unchanged. Evolve all four
with the same held elastic graph, cold seed and nominal DOP853
tolerance as the [first unsplit daughter controls](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-PREPARED-UNSPLIT-2026-11-27.md).

| Preparation | Mean cold source | Shift from native |
|---|---:|---:|
| Native | 0.06649893898 | 0 |
| Local daughter | 0.06701618429 | +5.17245309e-4 |
| Outside daughter | 0.06572727839 | -7.71660590e-4 |
| Entire daughter | 0.06624468123 | -2.54257750e-4 |

The positive local effect offsets about 67% of the outside reduction.
The four-arm interaction is +1.5753e-7 in the first-window mean, about
0.062% of the absolute whole-daughter shift. Thus the two input changes
are nearly additive at this finite held-graph comparison, but of
opposite sign. This interaction is not the interaction between the
two elastic channels, which was evaluated under a different set of
controls. The whole-daughter reduction here is only about one third
of the outside reduction. The earlier [98.69% estimate](INDEPENDENT-EARLY-POLE-PREPARED-OUTER-HARD-SWAPS-2026-11-27.md)
for the fraction retained by the outside change is specific to the
16-step active split and does not carry over to this simultaneous solve.
The cumulative interaction also changes sign during the window.
All four outputs pass the existing conservation gates; the minuscule
interaction is a model and numerical diagnostic, not a universal bound.

From `kinetics/`, run the simultaneous reference with
`--prepared-input` set to the accepted split whole-daughter output.
Run `early_pole_unsplit_daughter_additivity_compare.py` with four
`--raw-inputs` and four `--histories` in the order native, local,
outside, entire. Its input checks require the 2 and 18-cell
perturbations to be disjoint and exactly reproduce the whole input.
