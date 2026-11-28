# Cross time resolution and outer support on four completed corrected runs

Extend the [matched corrected ladder](INDEPENDENT-WARM-CORRECTED-TIME-AND-SUPPORT-2026-11-26.md)
with a q64 trajectory using 48 time intervals rather than 24.
It starts from the same q64 seed as the coarse widened run and
keeps the wider internal targets and its own accepted predecessor
at every step. The archived q32-48 states are only nonlinear
initial guesses; all 48 guided solves accept on q64 without
fallback. The two q32 and two q64 trajectories end at the same
7.394049 MeV after expansion by 1.10.

At every common clock point, let d32 and d64 be the signed
24→48 changes in the comoving number, full energy and pressure
per shared momentum shell. Compute d64−d32 before taking the
shell L1. This is the mixed time/support difference, not the
difference of two absolute errors. Normalize each metric to
its initial represented q32 moment by species.

| Metric, maximum over shared times and moments | Relative shell L1 |
| --- | ---: |
| Time refinement q32 | 5.8093e-4 |
| Time refinement q64 | 5.8093e-4 |
| Outer support q32→q64 at 24 intervals | 3.9739e-6 |
| Outer support q32→q64 at 48 intervals | 4.0090e-6 |
| Signed mixed time/support difference | 3.5612e-8 |

The q64-48 cold-cell average declines at all 48 intervals and
ends at 4.09445744e-5. Its last slope is −5.09616e-5 per ln(a), and
the explicitly computed first-cell prompt channel is still open.
The q32 counterpart has a geometrically closed first-cell prompt
channel. The absolute mixed effect on the cold occupation itself
is below 2.77e-11 at the common clock points. This comparison
separates the two selected discretization changes. It says
nothing about joint convergence outside these two represented
supports, or convergence of missing interactions.
Replaying that last finer q64 step gives the exact archived
hash and apportions the cold decrease per ln(a) into
−1.00895e-5 pair, −4.08669e-5 below-pole elastic, and −5.27337e-9
prompt. Elastic still accounts for about 80% locally.

The resolved q64 positive high inverse energy rate at the late
cut .70779 GeV is 4.1173e-32 of retained dark energy times H at
48 steps, compared with 4.1300e-32 at 24 (−0.31%). This is the
[frozen support diagnostic](INDEPENDENT-WARM-HIGH-INVERSE-FROZEN-PRESSURE-TEST-2026-11-26.md),
not an integrated exposure or a physical missing-tail bound.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# Prepare ladder/{seed32.json,q32-24,q32-48,wide} in the linked ladder.
python -W error warm_q64_multishell_guided_continuation_reference.py EOS.dat \
  results/warm-recoil-crossing-q64-terminal-2026-11-26.json.gz ladder/q32-48 ladder/wide48 --blocks 48
python -W error warm_multishell_corrected_support_audit_reference.py EOS.dat \
  ladder/seed32.json ladder/q32-48 ladder/wide48 48
python -W error warm_multishell_corrected_time_support_interaction_reference.py EOS.dat \
  ladder/seed32.json ladder/q32-24 ladder/q32-48 ladder/wide ladder/wide48
python -W error warm_corrected_archive_high_parent_recheck_reference.py EOS.dat \
  results/warm-corrected-anchor-signed-floor-2026-11-26.json ladder/wide48/q64/step-48.json
python -W error warm_q64_corrected_actual_cold_budget_reference.py EOS.dat \
  ladder/wide48/q64/step-47.json ladder/wide48/q64/step-48.json ladder/q32-48/step-48.json
~~~

[Signed mixed differences](results/warm-corrected-time-support-interaction-2026-11-26.json.gz);
[48-interval support audit](results/warm-multishell-corrected-support-audit48-2026-11-26.json.gz);
[q64 step 47](results/warm-multishell-corrected-q64-fine-stage-47-2026-11-26.json.gz);
[q64 step 48](results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz);
[fine high-band probe](results/warm-corrected-archive-high-parent-q64fine-2026-11-26.json).
[Fine final cold-cell split ledger](results/warm-corrected-fine-terminal-cold-budget-2026-11-26.json).
