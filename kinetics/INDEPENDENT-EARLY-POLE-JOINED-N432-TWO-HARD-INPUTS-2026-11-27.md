# Both hard-input histories pass the tighter 432-step refinement

Finish the matched 7.2e−5 fixed-scale five-channel comparison with the [same 5e−15 per-step elastic budget](INDEPENDENT-EARLY-POLE-JOINED-N432-NATIVE-RECOVERY-2026-11-27.md) in both fine runs. Each fine and coarse pair uses identical initial hard moments, cold bin counts, supports, interaction amplitudes and active unlike screening. Only the time grid and elastic acceptance budget differ between the *accepted* 216/432 runs. The tighter acceptance at the same 432-step width changes the stored integrated scalar history by at most 1.05e−17 in either input, used only as a solver sensitivity because the earlier 432 ledgers were rejected.

| Hard input | Gain at 216 | Gain at 432 | Relative gain change | Terminal energy error at 432 |
| --- | ---: | ---: | ---: | ---: |
| Archived steps | 7.845251615e−6 | 7.847737198e−6 | +0.03168% | +2.0872e−14 |
| Count-conserving continuous | 7.833774694e−6 | 7.838063723e−6 | +0.05475% | +1.5321e−14 |

Terminal energy and entropy pass the unchanged 3e−12 and 3e−11 limits in both arms. For the continuous input, the default 432 attempt failed at +2.18681e−11; the tight signed sum of elastic errors is +1.7097e−14. The two accepted time refinements change the first coarse-width averaged cold rate by −4.99% and +2.12%, respectively, despite changing the terminal instantaneous rates by only about +.014%. The early rate is less resolved than the accumulated gain.

At the finer width, the integrated continuous-input gain is .1233% lower than the native-input gain, compared to .1463% at the previous width. Its last mean rate is .1554% lower. At the matched terminal snapshot, parent and daughter relative occupation L1 differences in the aperture are .193% and .402%; on full common support they are .731% and .546%. Thus the small final count difference does not establish spectral memory loss, freeze-out, or a continuum-limit rate. This calculation has no cosmological expansion.

From `kinetics/`, write checkpoints to fresh files:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for p in steps cell-conservative; do
  python -W error early_pole_cold_screened_split_reference.py fine \
    --duration 7.2e-5 --count 432 --pair --same-soft-strength 1 --compton-strength 1 \
    --conversion-strength 1 --screening-update full --projection "$p" --snapshots \
    --elastic-energy-budget 5e-15 --checkpoint "$p-tight.ckpt.json" \
    --checkpoint-every 12 > "$p-tight.json"
  python -W error early_pole_joined_split_checkpoint_ledger_reference.py "$p-tight.ckpt.json"
  python -W error early_pole_joined_temporal_refinement_reference.py \
    results/early-pole-joined-fine-216-memory-2026-11-27.json.gz "$p-tight.json" --projection "$p"
done
python -W error early_pole_joined_history_memory_reference.py steps-tight.json cell-conservative-tight.json
```

For the numerical sensitivity only, reproduce the rejected checkpoints without `--elastic-energy-budget`; the original terminal failure is expected. The companion `early_pole_joined_long_refinement_bundle_reference.py` assembles the accepted trajectories, signed same-width check, and both time/reconstruction comparisons from those checkpoints and accepted output. Raw: `results/early-pole-joined-fine-432-time-refinement-2026-11-27.json.gz`.
