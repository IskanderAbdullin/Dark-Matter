# Live moments with a frozen exterior log interaction

Combine the [live logarithmic energy](INDEPENDENT-FINITE-MASS-QED-LIVE-LOG-FRW-2026-11-27.md)
with the [active/exterior partition](INDEPENDENT-FINITE-MASS-QED-LOG-EXTERIOR-2026-11-27.md).
Initialize the two live fermion populations by integrating the fine
and third archive occupations into the same comoving two-moment packet
cells. The photon window and the scalar are the reference fixture's
inputs. Outside the active charged window, hold the represented
occupations fixed above the hard floor q=0.25 and let the explicit
comoving-energy factors redshift. There is no exterior collision
flux and no parameter fit. Time and scale factor are prescribed.

The full quadratic has active–active, active–exterior, and
exterior–exterior terms. The new active gradient comes from the
middle term. Include its energy and its explicit fixed-comoving
scale derivative in the same budget; derive both from the mixed
quadrature used for the moment gradient. The purely exterior log
self-energy has no derivative with respect to active moments and is
subtracted along with its scale work throughout this ledger. It
would have to be restored when using an absolute total energy as
input to Friedmann evolution. The factorized exterior contribution
is held fixed between comparators, including its explicit redshift.

The two projections start from the same archive moment data in each
three-arm comparison: factorized only, log on active–active pairs,
and log on both active–active and active–exterior pairs. This
identifies the response to the cross term without changing the
collision clock, packet cache, or initial occupations. The cross
gradient enters both decay and pair remaps. The derivative checks
vary scale or moments separately, and order 14 is compared against
order 24 in the mixed log integration on the initial state.

The finite-mass log lowers the final scalar gain by 7.30–7.32
parts per million relative to the factorized runs. The cross term
alone accounts for 1.53 parts per million relative to the active-only
log control, roughly 21% of the combined correction in both archive
starts. The joint energy defect after subtracting scale work is
below 7.5e−16 in all six tests and entropy rises. The explicit
cross-energy first variation and trace agree with finite differences
within 1.3e−23 absolute; the mixed quadrature comparison changes the
initial energy, gradient or trace by at most 2.35e−7 relative.
These are small paired numerical differences, not a change in an
inferred dark-matter abundance.

This is a represented energy/entropy closure and source sensitivity
test with reference packet weights. It does not resolve the exterior
population's evolution, missing soft modes, full dispersion/widths,
or a physical cosmological rate.

Replay from kinetics with the prior archive snapshots:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_live_log_exterior_reference.py
```
[Raw receipt](results/finite-mass-qed-live-log-exterior-2026-11-27.json).
