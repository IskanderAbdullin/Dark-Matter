# Accepted-step restart preserves the frozen solver replay

The [frozen-chart accepted solver](early_pole_joined_simultaneous_chart_reference.py) is expensive enough that a lost environment can discard a long run. A completed-run state is not the same as a restart of the adaptive integrator mid-window. This reference checks a smaller accepted-step boundary without changing any production run.

For the tested version of SciPy DOP853, after an accepted integration step and after all dense queries in that step have been emitted, the next step depends on the current time, state, derivative, suggested absolute step, unchanged tolerances and endpoint. Keep the original sampled states, accepted times and the grouping of dense output calls. A replacement constructor first reuses the saved derivative, restores the step and RHS evaluation count, and proceeds without evaluating the original physical RHS a second time at the restart time. The serialized payload records the tested Python/NumPy/SciPy versions, exact stepping and dense-output source fingerprints, solver settings, and an integrity hash.

Retaining the sample array shape alone was insufficient for the frozen test. The first experiment preserved each mode but changed the aggregate cold count because transposing a continuous record changed NumPy's summation order. Reconstructing the same dense-output blocks and joining them with the same column layout removes that difference.

The synthetic adaptive fixture uses 63 accepted steps, 89 sampled times and 1127 RHS evaluations. Restarts at accepted steps 1, 3, 8, 29 and 45 reproduce the uninterrupted result bitwise, and the reconstructed dense matrix agrees with an independent SciPy integration in both values and strides. A modified step size under the old hash is rejected. Restarting the synthetic fixture in a separate Python invocation also gives the same final fingerprint. The actual two-panel H1 frozen both-elastic fixture reproduces its raw diagnostics bitwise in-process after accepted step 1 (347 physical RHS calls, duration 1.25e−10); its retained energy and entropy ledgers are unchanged.

This is a deliberately version-pinned reference, not a promise that every SciPy version stores compatible state. It does not support solver events or arbitrary external model side effects. A durable **production** restart must also persist model provenance, RHS-owned counters and diagnostics, and any external scratch state. The already running heavy branches are not altered by this audit.

From `kinetics/`, replay:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error accepted_dop853_checkpoint_reference.py
python -W error early_pole_accepted_step_resume_audit.py
# Then, as independent processes:
python -W error accepted_dop853_checkpoint_reference.py --write-checkpoint /tmp/pole-step.json
python -W error accepted_dop853_checkpoint_reference.py --restore-checkpoint /tmp/pole-step.json
```
The last two commands give final SHA-256 `020ff033e509e8b26ded0fac3a9f4c21432a61eb04e4d74c98e704dd80b7d913` for the recorded runtime. The [synthetic receipt](results/toy-accepted-dop853-resume-data-2026-11-28.json) and [frozen receipt](results/early-pole-adaptive-step-resume-data-2026-11-28.json) are separate to avoid conflating cross-process and physical tests.
