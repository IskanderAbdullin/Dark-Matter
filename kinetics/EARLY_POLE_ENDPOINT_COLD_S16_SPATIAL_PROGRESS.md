# Whole S16 cold band and the next adjacent spatial step

The S16 evolved cold-scalar chart has passed both native quadratures, with 158 complete panels per q mesh over [1e−6, 1.5e−5], bound to the archived x=.05 state. The maximum relative mesh discrepancy among the nine whole-band coefficients is 2.03e−13. Charge remains divided by two on the fine archived state; the recorded pole and scalar masses are chart-specific. Compared with S8, the integrated S16 forward number, inverse number, and forward excess are lower by 0.0205700%, 0.0241546%, and 0.1397482%, respectively.

With each chart's own native cold-cell number and mean excess fixed, the independently bounded minimum normalized TV on the same 17 display bins is:

| Adjacent levels | Forward number | Inverse number | Forward excess |
| --- | ---: | ---: | ---: |
| S2/S4 | 2.022556% | 2.020395% | 1.871110% |
| S4/S8 | 0.451486% | 0.450987% | 0.435061% |
| S8/S16 | 0.110803% | 0.110683% | 0.108277% |

The S8/S16 forward-number value is about one-quarter of S4/S8 on that common display. On the pair-specific 60-bin display the corresponding three S8/S16 results are 0.130642%, 0.129926%, and 0.122483%. The S8/S16 intervals on either display are at most 9.19e−10 wide in absolute TV for both moment restrictions. The old nominal solver-residual check does not pass every S16 case; we retain results verified instead by a separately computed fractional support lower bound and projected positive-profile upper witness. The two meshes agree on the 17-bin diagnostics to 1.59e−13 absolute. All numbers here refer to finite held panel mixtures, not a continuum extrapolation, physical density, or transfer relation.

Separately, when using shared reference moments at the single S8 optical mass, all 60 native common parents carry the pure S8 reference. At S16 reference weight 0.25, 0.5 or 0.75, S8 cannot carry 40 of the means; at weight 1 it cannot carry 41. S16 carries them all. The maximal common contiguous fallback retains 36 parents and the worst pure-S16 discrepancy on S8 relative to its requested mean is 0.5511%. The geometry-only split plan preserves all 60 display bins by replacing 41 S8 panels with 82 child panels per mesh (141 panels after assembly). It does not assign reaction coefficients to those children.

From `kinetics/` (both original panel pair files must be retained):
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
arc=results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
python -W error early_pole_endpoint_continuum_cold_band_receipt.py "$p"16-endpoint-continuum-cold-{1,2}.jsonl.gz --archive "$arc"
python -W error early_pole_endpoint_cold_common_grid_spatial_ladder.py --independent-bracket \
 --s2 "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz --s4 "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --s8 "$p"8-endpoint-continuum-cold-{1,2}.jsonl.gz --s16 "$p"16-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{2,4,8,16}-x05-state-$d.json.gz --archive "$arc"
for ch in forward_number inverse_number forward_excess; do
 python -W error early_pole_endpoint_cold_profile_exact_finite_distance.py --independent-bracket --levels 8 16 \
   --left "$p"8-endpoint-continuum-cold-{1,2}.jsonl.gz --right "$p"16-endpoint-continuum-cold-{1,2}.jsonl.gz \
   --states "$s"{8,16}-x05-state-$d.json.gz --archive "$arc" --observable "$ch"
done
python -W error early_pole_endpoint_cold_shared_reference_reintegration_plan.py \
 --left "$p"8-endpoint-continuum-cold-2.jsonl.gz --right "$p"16-endpoint-continuum-cold-2.jsonl.gz \
 --scan "results/early-pole-s8-s16-targeted-reintegration-native-geometry-$d.json"
```
