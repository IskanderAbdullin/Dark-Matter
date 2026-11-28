# The native-step long temporal refinement fails its terminal energy gate

The native-step fine trajectory completes all 432 macro steps at exposure 7.2e−5 with the same controls as the [216-step run](INDEPENDENT-EARLY-POLE-FINE-EXTENDED-MEMORY-2026-11-27.md). However, it **does not pass** the terminal energy ledger gate: the independently replayed defect is +2.26925e−11, above the unchanged limit 3e−12. The entropy ledger defect is −1.516e−12, within its separate gate. The provisional gain stored in the checkpoint is not accepted as a convergence result.

The sum of signed elastic substep energy defects is +2.27065e−11. Every individual elastic defect is below 2.914e−13; the first cumulative elastic sum to exceed the full-run gate is step 307. The difference between the terminal energy defect and the summed elastic defects is −1.40e−14. This locates the accumulation in the elastic acceptance budget for this finite scheme, rather than suggesting that a terminal tolerance should be loosened. It does not identify a physically missing contribution. The independently computed continuous-arm result is not included here.

Reproduce the stopped numerical check from `kinetics/` using a fresh writable checkpoint name. A terminal failure in the command is expected; the checkpoint contains the completed trajectory and independent controls, even when terminal output is rejected.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_cold_screened_split_reference.py fine \
  --duration 7.2e-5 --count 432 --same-soft-strength 1 --compton-strength 1 \
  --conversion-strength 1 --screening-update full --pair --projection steps --snapshots \
  --checkpoint native432.json --checkpoint-every 12 --progress-every 12
python -W error early_pole_joined_split_checkpoint_ledger_reference.py native432.json
```
Audit: `results/early-pole-joined-n432-ledger-rejection-2026-11-27.json`. A valid long time-refinement claim requires repeating the trajectory with adequately apportioned per-step accuracy and still respecting the terminal gate.
