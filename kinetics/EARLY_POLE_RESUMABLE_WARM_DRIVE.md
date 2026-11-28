# Resume long warm histories at exact internal step boundaries

The warm paired-clock driver now has opt-in per-arm progress and restart files. The JSON snapshot is atomically replaced only after a completed time step. It retains completed and target steps, reached x, occupation, all accumulated pressure-work trace values, extent caps, and the existing physical check records. A checksum detects accidental changes. A second signature binds the initial state, every clock knot, charged and scalar geometry, stepper identity, and the application's explicit provenance digest. Do not reuse a provenance digest for changed source or numerical inputs. Use a separate checkpoint prefix for each writer and stage.

This is a numerical restart facility, not a new chart selection or change to the collision method. The default path does not checkpoint. An actual short integration was interrupted after six of sixteen steps; the resumed output state and all physical trace records were identical to the uninterrupted result. Incorrect clock, corrupted data, and inconsistent trace lengths are rejected. A completed arm can be reused without performing any additional collision step.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest -q test_early_pole_resumable_drive
# First compute a SHA-256 digest of an immutable input/code manifest and assign it to PROVENANCE.
python -W error early_pole_endpoint_higher_warm_gauss_history.py \
 --warm-order 8 --x-end .0075 --steps 2048 --state /tmp/g8-warm-start.json \
 --checkpoint-prefix /tmp/g8-own-startup --checkpoint-provenance "$PROVENANCE" \
 --checkpoint-stride 128
```

Restart by repeating the identical invocation. Keep the individual files for the 2048 and 4096 arms until the final paired report and state have passed their ordinary acceptance checks. Do not treat a partial arm as physical evidence. This is not a cryptographic authentication mechanism for externally supplied data.
