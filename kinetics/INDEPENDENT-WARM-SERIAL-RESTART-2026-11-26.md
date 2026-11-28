# Restart a prescribed warm calculation at outer block boundaries

This extends the [portable exterior archive](INDEPENDENT-WARM-EXTERIOR-RESTART-2026-11-26.md).
Long calculations now have a serial reference driver which preserves
the planned log-expansion grid and commits the actual positive-shell
shape, known exterior cohorts, and accepted-block ledger after every
outer block. Its status is emitted only after each checkpoint is
written. An interrupted block is recomputed, not counted as complete.
The work directory locks while a process is writing and rejects a
restart with different seed, EOS, grid, solver flags, or tracked
main reference source files. This is not an alternative physical
equation, and does not reconstruct missing cohorts in an old shape.

A short, known-cohort 45 MeV test ran one block, exited, and resumed
through the second block. An attempted restart on an incompatible
grid was rejected without modifying the saved checkpoint. Reopening
the completed grid reproduced the checkpoint without recomputing a
block. Compared to the uninterrupted two-block run with otherwise
identical settings, weighted shell number/energy agreed to
2.4×10⁻¹⁵ and the total-count ledger to 8.9×10⁻¹⁵.
The resolved shells, visible clock and archived exterior quantities
were all checked. This does not cover a long physical trajectory
or a change of the phase-space grid.

## Exact alternative prompt derivative

The optional `--local-prompt-jacobian` switches to the existing
[parent-local mathematical derivative](INDEPENDENT-LOCAL-PROMPT-JACOBIAN-2026-11-26.md)
with the same forward collision operator and residual. At the true q32
12 MeV seed, the one-block expansion factors 1.0002 and 1.03 gave
weighted shell differences 1.8×10⁻¹⁵ and 2.4×10⁻¹⁴.
In separate sequential cloud measurements the optional derivative
used about 14–15% less wall time but about 20–21% greater peak RSS.
The control records use a fresh launcher per measured process.
This is a two-point performance observation, not a recommendation
to switch when memory is constrained, nor evidence that a full
q32 warm solve is faster or complete. The default is unchanged.

## Reproduction

~~~sh
export OPENBLAS_NUM_THREADS=1
S=results/warm-exterior-restart-producer-shape-2026-11-26.json
python -W error warm_coupled_serial_restart_reference.py EOS.dat "$S" state \
  --stop 1.0002 --steps 2 --max-blocks 1
cp state/checkpoint.json first-checkpoint.json
python -W error warm_coupled_serial_restart_reference.py EOS.dat "$S" state \
  --stop 1.0002 --steps 2
python -W error coupled_warm_pair_continuation_reference.py EOS.dat "$S" \
  --stop 1.0002 --steps 2 --track-exterior --save-snapshots unbroken-shape.json
python -W error warm_serial_restart_audit_reference.py "$S" first-checkpoint.json \
  state/checkpoint.json state/final-shape.json unbroken-shape.json

python -W error warm_q32_prompt_performance_reference.py EOS.dat \
  results/50MeV-from-empty-pair-quantum-q32-shape-2026-11-26.json q32-measurements
~~~
