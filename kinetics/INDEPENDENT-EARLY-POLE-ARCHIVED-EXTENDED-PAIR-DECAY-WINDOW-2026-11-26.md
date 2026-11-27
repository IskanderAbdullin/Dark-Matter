# Extending the evolved hard support changes the coupled scalar response

The [frozen archive budget](INDEPENDENT-EARLY-POLE-ARCHIVE-HARD-PAIR-EXACT-SUPPORT-2026-11-26.md)
shows that the first live hard box omits most of the daughter heat.
Here, use the same two archived x=.05 charged spectra to initialize
shared live photon/fermion cells at several upper momentum limits.
The archived count and centered energy of each live fermion cell
are integrated across the original edges and projected exactly
onto the two-moment closure. Exterior hard counts stay frozen and
their true positive moment remains in the reciprocal energy.
The difference between reconstructed and original initial field
moments is reported, rather than hiding it in a signed spectator.

Only a single scalar bin evolves in this comparison; the other scalar
bins have fixed counts but retain their energy and influence on the
live Hartree gap. The initial photon bath is unit-temperature BE.
Pair weights use the vacuum g⁴/H_ref normalization over the selected
domain at g=.005, with variational on-shell remap and triangle check.
An imposed scalar clock and a separate scalar vacuum-width-reference
clock both keep the same sign over the short reference exposure
8e−5. The latter reduces algebraically to the independent Yukawa
mode kernel in the no-medium Lorentz-invariant geometry. It still
omits medium amplitudes, widths, longer redshift and feedback into
exterior particles.

| Evolved upper limits photon / fermion | Fine relative suppression of scalar gain | Third relative suppression |
| --- | ---: | ---: |
| 2 / 2.2 | 3.77% | 3.65% |
| 4 / 5 | 7.94% | 7.70% |
| 8 / 10 | 8.55% | not run in this comparison |

The table uses the imposed scalar reference clock, and each row is
relative to its own decay-only control. Switching to the vacuum
scalar clock gives 3.77%/3.65% and 7.94%/7.70% in the first two
rows to displayed precision. The absolute gain of the chosen
native scalar bin, however, differs strongly across the archives,
whose bin edges and initial occupancies are not the same.
At the narrow fine support, incident refinement changes the
suppression from 3.77% to 3.80%. These are finite-box response
differences, not a convergence statement for the physical spectrum.

All accepted trajectories integrate energy, entropy, and total
photon-plus-charged number with separate ledgers. The instantaneous
collision checks are scaled by total absolute channel flux rather
than an absolute threshold calibrated for the much slower old toy
clock; the strict integrated energy and count checks are unchanged.
The widened model does not yet evolve the outgoing lower-than-.45
population. Revisiting that part and putting both scalar archives
on matched momentum support remain necessary.

Reproduce with early_pole_archived_extended_pair_decay_window_reference.py,
--grid fine/third, optionally --wide, --wider, --vacuum-clock or --level 1.
Raw: results/early-pole-archived-extended-pair-decay-window-2026-11-26.json.gz.
