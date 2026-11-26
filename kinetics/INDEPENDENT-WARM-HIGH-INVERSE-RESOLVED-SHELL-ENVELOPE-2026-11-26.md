# Bound the retained high inverse rate without trusting a q64 reconstruction

The [known-cohort birth-moment bound](INDEPENDENT-WARM-HIGH-INVERSE-KNOWN-COHORT-SHAPE-BOUND-2026-11-26.md)
concerned arbitrary shapes of stored exterior receipts at
activation. Here instead use only the number and energy in
each explicitly retained q64 shell, at the archived 9.748
and 9.027 MeV states. Hold their represented support fixed.

For each scalar–fermion shell pair, their two upper
energies must sum to the required minimum parent energy.
Each upper edge also imposes a necessary lower energy on
the other shell. Bound the eligible number with the stored
cell number and, when partial, the energy excess above the
lower shell edge. Divide by the minimum (pE) on the
eligible part to bound the radial energy-occupation measure.
Use the maximum possible parent energy and the exact
isotropic bare inverse pre-factor. Dropping the angular
condition and final heavy Pauli blocking can only increase
this positive source.

At 9.027 MeV the fine q64 two-moment bound is
(9.21	imes10^{-30}) in units of current retained energy
times Hubble, versus (2.68	imes10^{-32}) for the archived
entropy reconstruction. The number-only bound at that
epoch is (2.85	imes10^{-29}). The cut is 1 GeV in the
original 10.391 MeV comoving reference. The tested bounds
remain small without choosing a particular sub-shell
entropy profile. They are not bounds on inherited untracked
exterior modes, other heavy reactions, subsequent re-entry,
or the whole cosmological history.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_high_inverse_resolved_shell_envelope_reference.py EOS.dat
~~~

[Diagnostics](results/warm-high-inverse-resolved-shell-envelope-2026-11-26.json).
