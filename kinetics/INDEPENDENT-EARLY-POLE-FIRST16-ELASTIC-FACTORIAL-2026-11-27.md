# Finite early elastic effects cannot be added from their initial slopes

Use the matched first interval of 3.3333333333333335e−7 and the common 16-step grid. Keep reciprocal scalar decay, photon pairs, Compton and flavor conversion enabled. Toggle the screened unlike- and same-flavor terms separately, refreshing the unlike response only when that channel is present. Within each hard-input lineage the prepared moment state, cold counts, support and initial instantaneous source are identical in all four cases. All cases pass the original terminal energy and entropy gates.

| Elastic terms enabled | Native archived steps | Count-conserving continuous |
| --- | ---: | ---: |
| Neither | .07862183007 | .05689138986 |
| Unlike only | .06670942048 | .06207451735 |
| Same only | .06733209069 | .06022773304 |
| Both | .06694243427 | .06528878623 |

The signed interaction in the mean rate, defined as both minus unlike minus same plus neither, is +.01152275 for the native input and −.00012207 for the continuous input. In the native case each channel alone suppresses production over this window, but it would be wrong to add the two suppressions: the fully coupled mean is close to each one-channel mean. In the continuous case the same-flavor term alone **increases** the finite mean, even though its [starting directional derivative](INDEPENDENT-EARLY-POLE-INITIAL-COLD-RATE-DIRECTION-2026-11-27.md) is negative. Thus the local slopes explain how the prepared time-zero states begin moving, not the accumulated contribution of each elastic term. This is a finite numerical comparison, not a channel decomposition along the fully coupled trajectory or a late-time prediction.

### The continuous interaction reverses sign inside the window

Apply the same factorial to the cumulative count at each saved step. For the native input all 16 interval contributions are positive; only 10.38% of the final interaction has accumulated at the first quarter. For the continuous input the interval contribution becomes negative at step 8, and the cumulative value changes sign between 2.083e−7 and 2.292e−7. That reversal is not represented by the small negative full-window average.

In a smooth simultaneous flow, let R be the instantaneous scalar source, and U and V the two elastic vector fields. The common R and first source derivative cancel in the toggle difference, so the first possible integrated mixed term is cubic, with coefficient `[2 D²R(U,V) + DR(DU·V + DV·U)]/6`. By themselves the stored first two prefix ratios, 3.92 and 4.74, do not resolve the nonzero cubic asymptote (which would give a ratio of eight). Finite split width and changing support can matter, so no initial mixed coefficient is estimated from this grid.

From `kinetics/`, with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`:
```sh
for p in steps cell-conservative; do
  for cfg in none unlike same both; do
    b=(); u=0; s=0; screen=frozen
    case "$cfg" in
      unlike) u=1; screen=full ;;
      same) s=1 ;;
      both) u=1; s=1; screen=full ;;
    esac
    if test "$cfg" != none; then b=(--elastic-energy-budget 5e-15); fi
    python -W error early_pole_cold_screened_split_reference.py fine \
      --duration 3.3333333333333335e-7 --count 16 --pair --strength "$u" --same-soft-strength "$s" \
      --compton-strength 1 --conversion-strength 1 --screening-update "$screen" --projection "$p" \
      --snapshots "${b[@]}" > "$p-$cfg.json"
  done
  python -W error early_pole_initial_elastic_factorial_reference.py --projection "$p" \
    --none "$p-none.json" --unlike "$p-unlike.json" --same "$p-same.json" --both "$p-both.json"
done
```

Raw: `results/early-pole-first16-elastic-factorial-2026-11-27.json`.
