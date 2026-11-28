# Fixed projections separate raw band totals under shared scalar moments

On the 17 shared S2 in-band scalar parents, give both charts their same S2 count and mean excess at the same definition mass. Split a cell only among its archived q panels, preserve positivity, and leave its own optical/charge response unchanged. F, B and X are the forward-number, inverse-number, and forward-excess totals divided by their respective original S2 totals. Each aggregate range overlaps on its own. Those overlaps cannot be combined independently: a reweight changes all three totals together.

For 0.69030 F + 0.71998 B − 0.07156 X, the S2 range [1.33871438939648, 1.33872149508632] is disjoint from the S4 range [1.33860724546779, 1.33862321950580]. Therefore some band total must differ by at least 6.15248e-5 of its S2 reference total. Even the two number channels alone separate: for 0.90785 F − 0.41929 B, the gap is 1.17497e-5, forcing a difference of at least 8.85341e-6 in one reference-normalized number rate. Both projections overlap if only counts are prescribed; they then do not give any lower bound on agreement.

Once directions are fixed, neither check uses an optimizer. The scalar parents are independent; the minimum or maximum of a linear functional is the sum of the respective extrema over their stored convex vertices and the fixed spontaneous part. The lower bound on the maximum normalized difference divides the projection gap by the sum of absolute coefficients. Both accepted meshes agree to 1.5e-13 on the interval endpoints. This is a finite rate-amplitude diagnostic for archived states, not proof of continuous support separation or physical rates.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar; d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_cold_shared_aggregate_separation.py
python -W error early_pole_endpoint_cold_shared_aggregate_separation.py \
 --left "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --right "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{2,4}-x05-state-$d.json.gz \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
```
