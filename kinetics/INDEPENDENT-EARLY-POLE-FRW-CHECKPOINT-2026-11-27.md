# Checkpoint recovery for the joined expanding reference

The [joined screened FRW reference](INDEPENDENT-EARLY-POLE-JOINED-SCREENED-FRW-2026-11-26.md) now optionally saves a complete state at the boundary between macro steps. This includes the redshift-trace and entropy ledgers, the cold count history, and the implicit step diagnostics. The solver validates the numerical controls, initial state, discretization and stored state digest before resuming. Incomplete split steps are not written. Saving is opt-in and the same physical flags must be supplied on resume; this does not alter the collision or expansion algorithm.

A five-step expanding fine run with the conditional full H clock agrees byte-for-byte with an uninterrupted run after deliberately stopping following step two and resuming. A different supplied duration is rejected. Both continuation and the normal run retain the existing conservation and entropy checks.

From `kinetics/`, use a new checkpoint path:
```sh
export OPENBLAS_NUM_THREADS=1
args=(fine --duration 2e-6 --steps 5 --expands --full-H)
python -W error early_pole_joined_screened_frw_reference.py "${args[@]}" \
  --checkpoint resumable-frw.json --stop-after 2 > stopped.json
python -W error early_pole_joined_screened_frw_reference.py "${args[@]}" \
  --checkpoint resumable-frw.json --resume > resumed.json
python -W error early_pole_joined_screened_frw_reference.py "${args[@]}" > full.json
cmp resumed.json full.json
```
The reference final output in `results/early-pole-joined-frw-resume-check-2026-11-27.json` has SHA-256 `4cf03718fd79e6d564ee46ce3458bdf0322a93449d9f529506d2f061b4f08e14` on this tested setup.
