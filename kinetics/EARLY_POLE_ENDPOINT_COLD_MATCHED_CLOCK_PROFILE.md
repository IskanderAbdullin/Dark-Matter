# The S2–S4 cold shape discrepancy survives a matched final time refinement

Use each level's own independently evolved x=.04 restart, then compare explicitly paired 1024- and 2048-step paths to x=.05. Charge geometry stays at D2. Both q meshes accept the full band [1e-6,1.5e-5] at all four endpoints. The state hashes in the profile receipt link the independently gated integral receipt; the normalized comparison uses the exact intersection of all four edge sets (17 bins), with no interpolated responses.

| Quantity, held on each final clock | 1024 final steps | 2048 final steps |
| --- | ---: | ---: |
| S4/S2 integrated forward number response minus one | −0.331787% | −0.330178% |
| Unreweighted forward-number TV | 2.023986% | 2.024020% |
| Minimum forward TV keeping each chart's own counts and means | 2.022633% | 2.022556% |
| Minimum inverse TV keeping those moments | 2.020530% | 2.020395% |
| Minimum forward-excess TV keeping those moments | 1.873357% | 1.871110% |

Within a single chart, changing only the last time step shifts the direct forward amplitude by about −0.28% (coarse relative to fine). That shift largely cancels in the matched **ratio**. Allowing positive within-cell reweights at each chart's own retained moments still leaves the shape discrepancy. The fine-minus-coarse forward conditioned TV change is only −7.75e-7 in probability fraction. Mesh-to-mesh maximum reported difference 1.49e-13. This rules out the final step size as an explanation of the observed S2–S4 shape separation, but does not constrain earlier time-history error, higher scalar levels or physical/continuum convergence. The optional higher-level ladder interfaces do not imply any accepted higher-level bands.

From `kinetics/`, use the native and coarse q panel archives, the four associated states, and the existing independent optical reference (`results/finite-mass-qed-fixed-pole-cold-scalar-band-2026-11-28.json`). To reproduce both receipts:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
coarse2="$s"2-coarse-x05-state-$d.json.gz
coarse4=results/early-pole-geometric-warm-w1-coarse-x05-state-$d.json.gz
base=results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
args=(--s2-coarse "$p"2-coarse-clock-endpoint-continuum-cold-{1,2}.jsonl.gz
 --s2-fine "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz
 --s4-coarse "$p"4-coarse-clock-endpoint-continuum-cold-{1,2}.jsonl.gz
 --s4-fine "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz --archive "$base")
python -W error early_pole_endpoint_cold_scalar_band_clock_contrast.py "${args[@]}" \
 --states "$coarse2" "$s"2-x05-state-$d.json.gz "$coarse4" "$s"4-x05-state-$d.json.gz > /tmp/matched-cold-integral.json
python -W error early_pole_endpoint_cold_matched_clock_profile.py "${args[@]}" \
 --states "$coarse2" "$coarse4" "$s"2-x05-state-$d.json.gz "$s"4-x05-state-$d.json.gz \
 --integrated-clock-receipt /tmp/matched-cold-integral.json
```

The archived S2 coarse endpoint can also be regenerated with `early_pole_endpoint_s2_cold_clock_recovery.py` using the actual .04 state, .05 fine state, and .05 paired report, then an output argument. The default ladder still accepts its original S1/S2/S4 inputs.
