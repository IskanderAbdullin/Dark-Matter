# Restart a cold split run at a completed macro step

The split solver can optionally write a checkpoint using a temporary file followed by atomic replacement. This is an operational change, not a new term in the dynamics. The stored state contains the hard and scalar coordinates, integrated entropy/work and all output histories, including response updates, snapshots and rejection counters. Resume reconstructs the fixed geometry from the same inputs and starts at the next completed macro boundary. It does not change the interval, step count, numerical version or model parameters. No two jobs should write to the same checkpoint; omission of all checkpoint options preserves the usual output.

The two-step fixture enables photon pairs, hard Compton, hard flavor conversion, unlike screened elasticity and the same-flavor screened pole. Stopping after the first macro step and resuming produces the same serialized result hash as the uninterrupted run, including the final gain, ledger and occupation snapshots. Explicit invalid state content is rejected. Checkpoints are snapshots of a single run, not published science or a proof of convergence.

From `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
python early_pole_split_restart_fixture_reference.py
# For a long run, add these to the normal solver command:
# --checkpoint /path/to/unique-run.json --checkpoint-every 8 --progress-every 8
# Repeat the same physical/numerical arguments and add --resume after interruption.
```

Progress is written to stderr; completed JSON results remain on stdout. The optional `--stop-after` is for deterministic partial-run tests. The verified fixture record is `results/early-pole-split-restart-fixture-2026-11-27.json`.
