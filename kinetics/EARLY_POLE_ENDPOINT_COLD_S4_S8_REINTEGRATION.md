# Direct S4 reintegration restores all 31 S4/S8 conditioning cells

The archived geometry scan showed that 13 of the 31 common cold parents could not carry the requested S8 mean with nonnegative weights on the native S4 panel family. We now integrate both children of each of those 13 S4 panels from the same evolved S4 snapshot, in both angular/radial meshes. The unmatched old panels are carried through without reevaluation. No child coefficient is inferred from its parent and S8 is unchanged.

The joined 84-panel pairs pass the retained full-band gate. Comparing the sums of the children with the corresponding original S4 panels in nine components gives maximum relative differences of 2.81e-13 and 1.84e-13 on the two meshes. The whole-band paired component discrepancy is 1.64e-13. These gates test partition consistency at the fixed S4 state, not convergence of the evolved spatial grid.

All 31 parents support the count and mean in both charts for each of five reference mixtures. Relative to the old 22-parent fallback, the additional constraints raise the scale-normalized joint minimum. In particular, the pure S4 reference recovers the previous full-parent result, while the pure S8 reference is now feasible without fallback.

S8 weight | Shared count only, minimum L1 | Shared count and mean, minimum L1
---:|---:|---:
0.00 | 0.000046512206 | 0.000133596120
0.25 | 0.000046494675 | 0.000133658939
0.50 | 0.000046477144 | 0.000133721757
0.75 | 0.000046459613 | 0.000133784576
1.00 | 0.000046442086 | 0.000133847395

All loss entries use the same fixed S4 channel scales. The independently evaluated support/simplex bracket widths for count-and-mean minima in this table are at most 6.69e-13; the meshes agree on the reported optima to 2.35e-13. These are floating-point checks, not directed-rounding intervals. This comparison concerns finite nonnegative mixtures of held reaction panels, not unconstrained continuum spectra, a density prediction, or a transfer function.

## Reproduce from kinetics/

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28
p=results/fine-charge2-scalar
s=results/early-pole-fine-charge2-scalar
old="$p"4-endpoint-continuum-cold
right="$p"8-endpoint-continuum-cold
joined="$p"4-s8-targeted-reintegrated-endpoint-continuum-cold
plan=results/early-pole-s4-s8-targeted-reintegration-split-plan-$d.json
arc=results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
work=$(mktemp -d)
gzip -dc "$s"4-x05-state-$d.json.gz > "$work/state.json"
a=early_pole_endpoint_continuum_cold_band_pilot.py
for i in 1 2; do
 python -W error "$a" "$work/state.json" "$i" --topology-only \
   --mirror-azimuth --split-plan "$plan" > "$work/topology-$i.json"
 python -W error "$a" "$work/state.json" "$i" "$work/patch-$i.jsonl" \
   --mirror-azimuth --split-plan "$plan" --patch-only
 python -W error early_pole_endpoint_cold_targeted_patch_join.py \
   "$old-$i.jsonl.gz" "$work/patch-$i.jsonl" "$work/topology-$i.json" "$work/join-$i.jsonl"
 gzip -dc "$joined-$i.jsonl.gz" | cmp "$work/join-$i.jsonl" -
done
python -W error early_pole_endpoint_continuum_cold_band_receipt.py \
  "$joined"-{1,2}.jsonl.gz --archive "$arc"
python -W error early_pole_endpoint_cold_shared_reference_scan.py --levels 4 8 \
  --left "$joined"-{1,2}.jsonl.gz --right "$right"-{1,2}.jsonl.gz \
  --states "$s"{4,8}-x05-state-$d.json.gz --archive "$arc" \
  --left-baseline "$old"-{1,2}.jsonl.gz --left-split-plan "$plan"
```
