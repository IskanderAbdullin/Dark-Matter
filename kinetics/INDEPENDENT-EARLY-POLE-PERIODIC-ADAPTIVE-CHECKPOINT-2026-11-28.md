# Optional periodic solution checkpoints across Python processes

The [accepted-step reference](INDEPENDENT-EARLY-POLE-ACCEPTED-STEP-RESUME-2026-11-28.md) tested the solver state separately. The joined chart now has a strictly opt-in checkpoint mode. Its default integration entry point is unchanged. There is no automatic switch on a run that was already launched.

The opt-in mode writes at accepted steps only, after any report points for that step. The file includes the adaptive time, state, derivative, suggested next step, sample blocks and SciPy counters, along with the frozen RHS diagnostics: call count, reconstructed coordinate errors and sampled condition numbers. Adaptive source and numerical module fingerprints, initial archive identity, scope and the exact report grid are validated before restore. The sidecar is atomically replaced and synced. It is not silently accepted as the legacy completed-segment end state, which is a different operation and may shift the window.

A fresh-process H1 test with 32 modes, both elastic channels, duration 1.25e−10 and two panels reproduced the uninterrupted result byte for byte. The final file SHA-256 is `c3ec1f3426ef81114fd1d8a6fb8180bfe23dde4ac42a9a46ff10ddf6232e1bb7`. The last intermediate sidecar captured 27 accepted steps, two report samples and 332 physical RHS calls; the reconstructed result has the same 347 calls, the same per-mode values and aggregate counts, and unchanged energy and entropy diagnostics. The sidecar was written while the periodic run continued, then used by a separate new Python process. This verifies the exact scenario tested, not every nondefault process or an H8 two-window trajectory.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
seed=results/early-pole-first-window-hard1-cold32-paired-initializer-2026-11-27.json.gz
args=("$seed" --projection steps --duration 1.25e-10 --panels 2 --rtol 3e-9 --base-atol 2e-14 --record-mode-history)
prog=early_pole_joined_simultaneous_chart_reference.py
python -W error "$prog" "${args[@]}" > direct.json
python -W error "$prog" "${args[@]}" --adaptive-mode run --adaptive-path h1-step.json > periodic.json
python -W error "$prog" "${args[@]}" --adaptive-mode resume --adaptive-path h1-step.json > resumed.json
cmp direct.json periodic.json
cmp periodic.json resumed.json
```

Use a different sidecar per branch and per intended solver scope. To guard against loss of an execution environment, retain validated sidecars somewhere durable before the workspace is reset; an atomic write alone cannot survive a lost volume. The [physical checkpoint receipt](results/early-pole-h1-periodic-process-resume-2026-11-28.json) and [synthetic periodic check](results/toy-accepted-dop853-resume-data-2026-11-28.json) do not claim that any existing H8 run was checkpointed or rerun.
