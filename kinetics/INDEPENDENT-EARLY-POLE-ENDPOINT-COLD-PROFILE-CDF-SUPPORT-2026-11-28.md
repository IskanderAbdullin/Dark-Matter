# What the cold profiles do and do not imply about momentum displacement

A positive common-bin total variation minimum does not translate directly into a large displacement in q: edges of neighboring panels touch. Use the same fully resolved S2/S4 panel amounts, with optional finite positive reweights at fixed cell count and mean excess. Unlike TV accounting on shared bins, here deliberately allow **arbitrary within-panel positions of the integrated response**, without recovering its detailed shape.

At each open interval on the union of exact panel edges, all completely preceding panels must have contributed to the CDF. Any intersected panel can contribute anything between zero and its whole amount. Minimize the first cumulative normalized weight and maximize the second using the correlated finite moment reweighting for each chart. Separate the two outward CDF intervals only where nonoverlap is guaranteed throughout the open interval. Integrating the lower/upper CDF differences produces a conservative interval for the one-dimensional Wasserstein distance (W_1) in q units.

| Allowed normalized profile | Lower (W_1) bound | Upper bound |
| --- | ---: | ---: |
| Unchanged raw panel amounts | 4.155e−10 | 4.246e−7 |
| Finite reweights, counts only | 0 | 6.903e−7 |
| Finite reweights, counts and mean excess | 4.968e−11 | 4.575e−7 |

The count-only lower bound is zero, which does **not** exhibit a common distribution and does not override the 0.8833% finite minimum common-bin TV. At fixed counts and means, nine q intervals still have guaranteed CDF separation, but the certified momentum distance is small. These outward bounds allow different reweights at different cuts when forming envelopes, so are not claimed tight. The attempted alternative full interval coupling LP did not pass an independent marginal-closure check on these small numbers; no optimum from that method is used here.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_cold_profile_cdf_support
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
python -W error early_pole_endpoint_cold_profile_cdf_support.py \
 --left "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --right "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{2,4}-x05-state-$d.json.gz --levels 2 4 \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
```
