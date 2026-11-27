# Saved held-chart cold modes separate local rebound from outside suppression

Record the eight cold counts at each saved time for native,
two-cell local and compensation-only held simultaneous outputs.
All three scalar histories replay the previously accepted scalar
results bit for bit; summed individual cold counts agree with
their respective stored scalar histories. This is the original
cold grid, hard support, response and numerical tolerance.

The terminal local-minus-native difference is positive in all
eight modes; the terminal compensation-only difference is negative
in all eight. The local paired count first goes nonnegative
between 6.25e-8 and 8.3333e-8 in the highest-momentum bin, but
not until 1.25e-7 to 1.4583e-7 in each of the lowest five bins.
The two intervening bins cross between those saved bands. No
compensation-only mode shows a saved negative-to-nonnegative
transition. Exposures are in model units, not wall time.

All eight terminal local contrasts are larger in the simultaneous
chart than in the held 64-step split. The fifth and sixth bins
dominate their difference (about (6.71\times10^{-5}) and
(3.75\times10^{-5}) in mean gain). These comparisons measure
finite responses in the declared discretization, not a starting
temperature or a relic spectrum.

Use `--record-mode-history` with each original prepared case and
the same fine bundle, then compare the mode histories, prior
scalar histories and accepted split mode receipt using
`early_pole_unsplit_prepared_mode_compare.py`.
