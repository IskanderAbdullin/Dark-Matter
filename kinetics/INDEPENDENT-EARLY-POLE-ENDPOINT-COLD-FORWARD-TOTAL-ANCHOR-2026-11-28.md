# Exact forward-total matching can mask other cold response errors

On the same original archived S2/S4 finite panels with held, chart-specific optics, prescribe identical S2 in-band scalar counts and mean excess on both. The forward band total is not set to an arbitrary number: the admissible difference interval spans zero, about [−5.94e-6,+1.666e-4] in original S2 units. Require the scalar reweights to yield exactly equal forward totals, without separately normalizing either chart or refitting the optics.

| Objective at shared count and mean | Without forward-total equality | With forward-total equality |
| --- | ---: | ---: |
| Forward number profile L1 | 0.000033209085 | 0.000066418170 |
| Joint forward, inverse and forward excess L1 | 0.000147575396 | 0.001856101889 |

All channels are scaled by their own fixed original S2 totals. The price of exact forward calibration is a change in coupled constraints: in the joint optimum the total of the two number-profile discrepancies is about 0.000154, while the forward kinetic-excess discrepancy alone is 0.00170244. With shared **counts alone** both forward profiles can still match at numerical precision, but the joint loss rises from 0.0000237792 to 0.0000334869 when totals match.

As an independent numerical check, normalize nonnegative weights separately in each cell, then restore the forward-total equality by a convex blend with the opposite forward-extreme vertex pair. A support functional with bounded sign weights plus a multiplier on the exact total gives a matching lower objective. The bracket and two panel meshes agree within 5.25e-12 absolute. This is an error tradeoff for an explicitly constrained finite counterfactual, not the physical size of a collision-process error.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar; d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_cold_forward_total_anchor.py
python -W error early_pole_endpoint_cold_forward_total_anchor.py \
 --left "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --right "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{2,4}-x05-state-$d.json.gz \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
```

## Repeat after targeted reintegration

The published five-panel S2 reintegration makes every reference blend feasible on all 17 original cells; the response bins remain unchanged. With 0, 25, 50, 75, and 100 percent S4 projected reference weight, the forward band total difference always still spans zero. The equal-forward joint objective ranges from 0.00185610 at the S2 reference to 0.00181550 at S4, while the relaxed minimum ranges from 0.00014758 to 0.00014814. Agreement to 4.23e-12 across archived meshes. The penalty for exact scalar forward calibration is therefore not specific to choosing S2's original moments.

To run this version, set `--left` to the two reintegrated S2 band files and add `--left-baseline` with the two native S2 band files and `--left-split-plan results/early-pole-s2-targeted-reintegration-split-plan-2026-11-28.json`. All other inputs remain the same.
