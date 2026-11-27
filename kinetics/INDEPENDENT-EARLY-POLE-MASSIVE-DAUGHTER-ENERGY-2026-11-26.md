# Finite fermion-daughter mass in the conservative early calculation

The [energy-conserving live parent stencil](INDEPENDENT-EARLY-POLE-EXACT-ENERGY-PARENT-STENCIL-2026-11-26.md)
is generalized to a fermion daughter of fixed vacuum mass relative to
the parent, not a fixed value of its mass divided by temperature.
This changes the aperture through the exact two-mass discriminant,
the vacuum-matched scalar Yukawa factor, the fermion shell mean energy,
and the fermion contribution to the massive trace. The shell-averaged
energy and trace use the same phase measure. The two daughter changes
and fractional parent changes still balance projected energy and number.

Starting at x=0.005 with thermal parents and empty daughters, compare
at x=0.05 with self-consistent H and no momentum redistribution:

| Fermion mass / parent mass | Coarse yield | Refined yield |
|---|---:|---:|
| 0 | 0.000364269 | 0.000346891 |
| 0.06 | 0.000383539 | 0.000370191 |
| 0.30 | 0.000385152 | 0.000381765 |

The refined yields are about 6.7% and 10.1% above the massless null,
not a general statement that more daughter mass increases production.
Varying the mass changes several quantities together. At mass ratio
0.06, adding isolated parent mixing at rate 100 per radiation reference
expansion time gives 0.000425810. All the listed runs close total
continuity to within 3 × 10⁻¹¹ of the final daughter energy. The exact
zero-ratio rerun matches the previous massless refined solve to a
relative 2.3 × 10⁻⁸.

On a third matched grid with 97 scalar, 156 daughter, 176 parent
cells, mass ratios 0, 0.06, and 0.30 yield 0.000346861,
0.000364093, and 0.000380008. The relative enhancements become
+5.0% and +9.6%. The third grid changes the 0.06 result by −1.65%
and the 0.30 result by −0.46% relative to the previous grid.

At x=0.02 and scalar pole/T = 0.0005 an independent finite-pole
continuum rate test gives maximum active-row seed errors of 1.26%,
1.17%, and 0.88% for the three masses on the refined grid. The
phase-weighted seed errors are below 0.26%. These seed tests do not
bound the full many-mode stimulated evolution. The pointwise static
detailed balance test is discrete; Pauli bounds and exact projected
energy do not eliminate aperture quadrature or mean-field uncertainty.

~~~sh
python -W error early_pole_massive_daughter_energy_stencil_reference.py --ratio .06 --scalar 65 --fermion 90 --parent 72
python -W error early_pole_massive_daughter_energy_stencil_reference.py --ratio .30 --scalar 65 --fermion 90 --parent 72
python -W error early_pole_massive_daughter_energy_stencil_reference.py --ratio .06 --seed-only --scalar 65 --fermion 90 --parent 72
~~~

[Raw models and check states](results/early-pole-massive-daughter-energy-2026-11-26.json).
