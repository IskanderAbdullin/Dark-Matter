# Frozen and refreshed cold quadrature on already-evolved hard states

Use the accepted native and local-daughter held 64-leaf
histories. Replay the instantaneous cold source exactly at nine
stored times, starting from the original eight scalar bins.
Without evolving any state, treat each saved scalar count as
uniform in its own shell volume and integrate the source within
those eight shells. This assumes a different within-bin shape
than direct integration over the finer archived initial scalar
density; do not compare the two as if they were the same data.
Keep the saved total and the separately solved scalar pole,
the reconstructed charged occupations, and the same dressed
energy functional.

Locate the seven sorted angular cuts separately in each arm
and time sample. Their maximum motion over the saved times
is 7.56e-10 (native) and 7.38e-10 (local) in scalar q/T.
A second evaluation holds the union of both arms' initial
cuts. At 16 Gauss points per interval this startup-cut
evaluation differs from refreshed cuts by at most 2.64e-10
in the paired instantaneous rate. With refreshed cuts,
doubling the Gauss order from 8 to 16 changes individual
mode rates by at most 7.76e-13.

The refreshed volume result differs from the original
eight-volume-centroid paired source by up to 9.16e-6, at
startup. It preserves every sampled paired rate sign;
the instantaneous paired crossing remains bracketed by
[2.0833e-8, 4.1667e-8] in exposure. This is an instantaneous
source crossing, distinct from the later crossing of the
cumulative number contrast. No evolved response or time
integral is recomputed in this comparison.

From `kinetics/`, reconstruct the two held split paths using
`early_pole_prepared_hard_swap_trajectory_reference.py` with
`--hold-screening --steps 64` for native and local-daughter.
Pass the original archive and the two raw paths to
`early_pole_daughter_frozen_snapshot_cold_volume.py` at
`--order 8/16` separately, with and without
`--hold-startup-cuts`. Compare the four receipts using
`early_pole_daughter_frozen_snapshot_cold_volume_compare.py`.
