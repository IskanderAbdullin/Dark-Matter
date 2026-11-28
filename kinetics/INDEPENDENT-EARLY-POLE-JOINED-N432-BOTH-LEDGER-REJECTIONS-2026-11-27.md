# The independent continuous-input trajectory also fails its long energy gate

The [negative native-step result](INDEPENDENT-EARLY-POLE-JOINED-N432-LEDGER-REJECTION-2026-11-27.md) is not peculiar to that initial reconstruction. Keep the same five active channels, screening rule, fine hard mesh and duration 7.2e−5, changing only the parent/daughter input to the count-conserving continuous reconstruction. After 432 steps its energy defect is +2.18681e−11 versus the same 3e−12 gate. The independently summed elastic defects contribute +2.18630e−11 and first exceed the terminal gate at step 312 (native steps: 307). The entropy defect is −1.477e−12, within its gate. Both original long refinements are rejected for quantitative convergence.

Both matched reruns with stricter per-step acceptance are separate from these failed controls. The terminal energy gate is unchanged.

Reproduce the audit by following the [native command](INDEPENDENT-EARLY-POLE-JOINED-N432-LEDGER-REJECTION-2026-11-27.md) with `--projection cell-conservative` and a new checkpoint path. Raw: `results/early-pole-joined-n432-continuous-ledger-rejection-2026-11-27.json`.
