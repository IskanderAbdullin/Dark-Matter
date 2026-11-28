# Direct reintegration restores all 60 shared S8/S16 conditioning cells

The native S8 chart could carry the requested S16 mean on only 19 of the 60 common parents at the pure S16 reference. The maximal native contiguous fallback retained 36. At the same saved S8 state and optical fit, replace each of the 41 specified original S8 panels with two directly evaluated children. S16 does not change. This takes 82 new coefficients on each of two angular/radial meshes and yields 141 complete S8 panels. No new coefficient is assigned from a parent integral.

The largest original panel recomposition defects over the nine channel/moment entries are 3.90e−13 and 2.58e−13 on the two meshes. The whole-band paired discrepancy is 6.25e−14. Count and mean can now both be held fixed in all 60 common parents, across the two charts and all five reference mixtures. The jointly constrained minimum three-field L1 uses the fixed original S8 scales:

| S16 fraction in moment reference | Count only | Count and mean |
| ---: | ---: | ---: |
| 0.00 | 3.855500e−6 | 8.797013e−6 |
| 0.25 | 3.855309e−6 | 8.771919e−6 |
| 0.50 | 3.855118e−6 | 8.746826e−6 |
| 0.75 | 3.854927e−6 | 8.721732e−6 |
| 1.00 | 3.854735e−6 | 8.696638e−6 |

The independent support/simplex bracket widths for the shown count-and-mean optima are below 3.78e−14; the paired mesh optimum differences are below 3.25e−13. The 36-parent native fallback objective and the restored 60-parent objective are different conditioned problems. The split evaluates the same fixed S8 snapshot, not a new evolved state and not a physical or asymptotic convergence test. It also does not align every actual S16 population edge; the latter requires a separate rate reintegration before both scalar-factor replacement orders can be exact.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
old="$p"8-endpoint-continuum-cold
right="$p"16-endpoint-continuum-cold
joined="$p"8-s16-targeted-reintegrated-endpoint-continuum-cold
plan=results/early-pole-s8-s16-targeted-reintegration-split-plan-$d.json
arc=results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
dir=$(mktemp -d)
gzip -dc "$s"8-x05-state-$d.json.gz > "$dir/state.json"
for n in 1 2; do
 python -W error early_pole_endpoint_continuum_cold_band_pilot.py "$dir/state.json" "$n" \
   --topology-only --mirror-azimuth --split-plan "$plan" > "$dir/topology-$n.json"
 python -W error early_pole_endpoint_continuum_cold_band_pilot.py "$dir/state.json" "$n" "$dir/patch-$n.jsonl" \
   --mirror-azimuth --split-plan "$plan" --patch-only
 python -W error early_pole_endpoint_cold_targeted_patch_join.py \
   "$old-$n.jsonl.gz" "$dir/patch-$n.jsonl" "$dir/topology-$n.json" "$dir/join-$n.jsonl"
 gzip -dc "$joined-$n.jsonl.gz" | cmp "$dir/join-$n.jsonl" -
done
python -W error early_pole_endpoint_continuum_cold_band_receipt.py "$joined"-{1,2}.jsonl.gz --archive "$arc"
python -W error early_pole_endpoint_cold_shared_reference_scan.py --levels 8 16 \
 --left "$joined"-{1,2}.jsonl.gz --right "$right"-{1,2}.jsonl.gz \
 --states "$s"{8,16}-x05-state-$d.json.gz --archive "$arc" \
 --left-baseline "$old"-{1,2}.jsonl.gz --left-split-plan "$plan"
```
