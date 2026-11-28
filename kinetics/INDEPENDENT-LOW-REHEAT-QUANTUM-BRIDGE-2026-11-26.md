# Carry the low-reheat quantum state into a resolved later spectrum

Restart from the reproduced 24-step low-reheat state created with
signed visible FD reactions, prompt return, below-pole elasticity
and stable scalar-pair reactions. Keep the actual shape and evolve
to a/a_archive=10 (about 0.3315 MeV photon temperature) with the
published pair-only radiation continuation. Preserve the earlier
exterior count as a collisionless contribution to the inclusive
comparison; it is not fed into the shell reactions.

At matched resolution the q16 and q128 continuations differ by
only 3.4e-10 in relative **inclusive** scalar count, while q128
retains virtually all newly produced late daughters. Extending
the grid changes the represented scalar streaming diagnostic
because particles at the upper end are included. No alteration
of the inner initial state is made.

The chained fixed-count abundance proxies including early exterior
counts are Ωχh²=0.11970159, Ωφh²=2.9911e-10. Late pair production
adds about 22.1% to the scalar count and removes about 0.0105%
of the initial inclusive stable count. The separate vacuum-source
48-step early bridge gave Ωφh²=3.0031e-10; the different time meshes
do not isolate a quantum correction. Both late segments omit
visible interactions, on-shell transfers and elastic exchange.

Using the emitted retained shell distribution, the mean collisionless
distance **from this later archive** to equality is 0.2030 kpc for
the stable fermion and 20.33 Mpc for the scalar on the q128 grid.
The q16 scalar distance is 19.80 Mpc. Kinematic RMS speed at equality
is about 1.295e-7 c for the stable fermion and 0.05084 c for the
represented wide scalar. These numbers extend the earlier
[geodesic comparison](INDEPENDENT-ARCHIVED-COLLISIONLESS-STREAMING-2026-11-26.md);
they are not the cumulative distance for particles produced
earlier, a transfer function or an abundance constraint.

## Reproduce the later stage

The [early command](INDEPENDENT-SCALAR-PAIR-FROM-EMPTY-2026-11-26.md)
can save the retained quantum shape with `--save-snapshots`.
From `kinetics` run:

```sh
export OPENBLAS_NUM_THREADS=1
python -W error scalar_pair_expanding_continuation_reference.py EOS.late.dat early_quantum_shape.json --stop 10 --steps 24 --pad 128 --save-snapshots late_wide_shape.json > late_wide.json
python -W error archived_collisionless_streaming_reference.py EOS.late.dat late_wide_shape.json
python -W error early_scalar_pair_budget_reference.py EOS.late.dat results/low-reheat-matched-vacuum-2026-11-26.json results/low-reheat-early-pair-vacuum-2026-11-26.json early_quantum.json --late late_wide.json --bridge-input-index 2
```
