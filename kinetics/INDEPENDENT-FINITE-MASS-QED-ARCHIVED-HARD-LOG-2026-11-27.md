# Finite-mass exchange on two frozen archive grids

Apply the [hard-shell finite-mass exchange kernel](INDEPENDENT-FINITE-MASS-QED-SHELL-LOG-2026-11-27.md)
to the actual represented endpoint parent and daughter occupations
in the previously published fine and third x=0.05 archived kinetic
snapshots. Treat each distribution as isotropic and charge symmetric,
with the archived bin constants retained. This is an on-shell real
correction on frozen inputs, not a rerun of their kinetic evolution.

Cross both archives with internal lower momentum cuts 0.1, 0.25,
and 0.5, cutting cells without rebinning or extrapolation. At external
momenta p=0.25, 0.5, 1, 2, 5 compute the same-momentum parent-minus-
daughter mass-squared shift. Integrate its diagonal logarithm with
weighted quadrature split at p. Separately compute the symmetric double
integral and its fixed-comoving trace derivative with shell quadrature
orders 20 and 36. The factorizable trace uses an explicitly assumed
thermal photon bath with the same momentum cut; the photon term
cancels out of the species difference at common p. We do not infer
a photon distribution from the archives.

* In all six grid/cut combinations and at all five tested external
  momenta the parent-minus-daughter shift stays positive; the log
  reduces it by about 0.013% to 0.379%. At the central cut 0.25 the
  reduction is 0.283–0.284% at p=0.25 and about 0.138% at p=1.
* The ratio of logarithmic to factorized hard interaction trace is
  12.4% in both grids at the central cut. It ranges from 8.2% to
  19.1% as the cut changes. This spread describes cut sensitivity of
  the restricted trace, not a bound on the full interaction trace.
* The two grid values of the factorized species difference still
  vary by about 0.57% at the central cut, and varying the lower
  momentum cut has a larger effect there than the log correction.
  A pointwise same-momentum comparison is not the change in occupied
  decay phase space, where the two fermions generally have unequal
  momenta. No conclusion on a threshold crossing follows.

Reproduce with the archived
[fine](results/early-pole-spectrum-fine-m06-2026-11-26.json.gz) and
[third](results/early-pole-spectrum-third-m06-2026-11-26.json.gz)
snapshots in the kinetics results directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_archived_hard_log_reference.py
```
[Raw receipt](results/finite-mass-qed-archived-hard-log-2026-11-27.json).
No thermal widths, screening resummation, matched decay collision
operator, soft extrapolation or cosmological inference is included.
