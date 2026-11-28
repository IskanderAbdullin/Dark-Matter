# Same-flavor elasticity dominates the early sign change in the prepared witness

Take the [same fine 97% parent witness](INDEPENDENT-EARLY-POLE-NATIVE-WITNESS-JOINED-RELAXATION-2026-11-27.md) and hold the 16 scalar bins, three hard aperture subdivisions and three other hard channels fixed. Compare the two screened elastic sectors individually against no elastic and the already joined result, following the same exposure and macro grid. Where present, the unlike response is refreshed independently along that branch; same-flavor screening stays held in all cases.

| Active elastic sectors | First sampled rate at 2.5e−9 | Final rate at 2e−8 | Accumulated scalar gain |
| --- | ---: | ---: | ---: |
| None | negative throughout | −.006817 | −1.36396e−10 |
| Unlike only | −.006283 | +.007238 | −1.18730e−11 |
| Same-flavor only | +.033400 | +.070312 | +1.10872e−9 |
| Both | +.030263 | +.060679 | +9.69149e−10 |

All four branches start at −.00682279. The small negative accumulated unlike-only result is near cancellation, even though its last two stored instantaneous rates are positive. The both-channel gain falls 2.64090e−10 below the sum of the two incremental gains over the no-elastic control, so adding the two separate changes overcounts the joined result.

An independent replay at the same times gives separate parent and daughter occupation integrals (here rates with the stimulated scalar factor). In the first joined macro interval the parent contribution rises from .026466 to .067948, while the daughter inverse magnitude rises from .033289 to .037685. In the same-only branch these become .067040 and .033640; in the unlike-only branch they become .031275 and .037558. Thus the fast net sign change in this example is mostly the replenished **parent** aperture in the same-flavor sector, partially opposed by growth of the daughter term. This is a measured finite-model attribution, not a statement about complete same-flavor rates.

An independent joined control doubles the number of hard aperture subdivisions to six, holding all other settings. Its initial rate is more negative (−.011013), but the first sample at 2.5e−9 is still positive (+.019227), and the net gain over 2e−8 is positive (+8.58286e−10). That differs by 11.4% from the three-subdivision result. The sign and occurrence of fast relaxation survive this refinement; a precise crossing time or positive gain magnitude is not hard-grid converged. Energy and entropy ledgers close on all finite tested runs. An arbitrary prepared perturbation does not establish how such profiles are populated.

From `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1
python early_pole_live_aperture_witness_trajectory_reference.py fine --subdivision 3 \
  --scalar-bins 16 --duration 2e-8 --mix .97 --pair --hard-slow
# Run with (--strength, --same-soft-strength)=(1,0), (0,1), (1,1).
python early_pole_cold_screened_split_reference.py fine --duration 2e-8 --count 8 \
  --strength 1 --same-soft-strength 0 --compton-strength 1 --conversion-strength 1 \
  --screening-update full --pair --snapshots --aperture-subdivision 3 --scalar-bins 16 \
  --witness-mix .97 --response-radial-order 32 > unlike.json
python early_pole_joined_snapshot_aperture_reference.py unlike.json
# Repeat (1,1) with --aperture-subdivision 6 for the higher hard-resolution control.
```

Raw: `results/early-pole-witness-elastic-sector-attribution-2026-11-27.json.gz`.
