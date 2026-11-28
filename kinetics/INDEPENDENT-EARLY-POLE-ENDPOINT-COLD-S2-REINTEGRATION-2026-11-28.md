# Targeted S2 panel reintegration recovers all 17 shared moment parents

Reintegrate both halves of five existing S2 q panels, selected solely because the projected S4 scalar mean was outside the attainable interval on its parent. The native S2 occupation, charge state, fitted optical parameters, whole q support and shared response bin edges do not change. The old, unaffected panel rows are reused exactly. The ten new evaluations per mesh use the same numerical orders and mirror reduction; each pair re-sums to its original nine coefficients within 1.66e-13. The full 62-panel quadratures agree to 1.34e-13. Both the split raw patches and their assembled outputs are archived. No panel response is fabricated by dividing a parent integral.

All 17 original moment parents now admit projected moments from both charts. At S4 reference weights 0%, 25%, 50%, 75%, 100%, the finite joint absolute L1 minima in common original S2 scale units are 0.000147575396, 0.000147608030, 0.000147718554, 0.000147931284 and 0.000148144612. The shared-profile obstruction survives after removing the moment realizability limitation, without dropping four conditioning boundaries as in the earlier 13-parent experiment. Each chart still retains its own response and finite panel occupation families. This is neither a physical rate bound nor a test of a newly evolved S2 trajectory.

From `kinetics/`, to replay the archived fresh patch assembly and compare the responses:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar; d=2026-11-28
r=results/early-pole-s2-targeted-reintegration-
base=results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
plan="$r"split-plan-$d.json
gzip -dc "$s"2-x05-state-$d.json.gz > /tmp/cold-s2-state.json
for j in 1 2; do
 python -W error early_pole_endpoint_continuum_cold_band_pilot.py /tmp/cold-s2-state.json "$j" --topology-only --mirror-azimuth --split-plan "$plan" > "/tmp/cold-s2-topology-$j.json"
 # Recompute the ten patch pieces if needed: same CLI, replace --topology-only
 # with --patch-only and add a fresh output JSONL positional argument.
 python -W error early_pole_endpoint_cold_targeted_patch_join.py \
   "$p"2-endpoint-continuum-cold-$j.jsonl.gz "$r"raw-patch-mesh-$j-$d.jsonl.gz \
   "/tmp/cold-s2-topology-$j.json" "/tmp/cold-s2-replayed-$j.jsonl"
done
python -W error early_pole_endpoint_cold_shared_reference_scan.py \
 --left /tmp/cold-s2-replayed-{1,2}.jsonl --left-baseline "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz --left-split-plan "$plan" \
 --right "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{2,4}-x05-state-$d.json.gz --archive "$base"
```
