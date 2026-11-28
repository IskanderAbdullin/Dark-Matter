# Sensitivity to the choice of shared reference scalar moments

Project both evolved, step-reconstructed scalar states onto the same 17 S2 in-band conditioning parents. Integrate each exact overlap in q² volume at the same reference S2 mass. Mix the two integrated count-excess pairs, not the means themselves. The original rate evaluation remains fixed: chart-native optical fits and charge response, unchanged accepted q meshes and comparison response bins.

When both charts get the same counts only, shifting the reference continuously from S2 to S4 changes the three-field joint finite optimum from 2.3779e-5 to 2.4369e-5 over the sampled five settings. Adding mean excess on **all** 17 parents changes what the finite support can realize. With S2 moments there is a valid solution. With any positive S4 fraction, the S2 panel in [1.09496e-6,1.31638e-6] can no longer match its moment: it has only one panel mean, while the weighted reference mean moves strictly above it. At 25% S4 there are four infeasible S2 parents, at 50% through 100% five, while S4 remains feasible. At 100% the largest relative mean miss on S2 is about 4.67%. This is an obstruction of the saved finite panel basis, not evidence for a continuum lower bound.

A dynamic program finds the greatest number of retained conditioning parents for which both projected endpoint moment pairs are feasible on both finite panel sets. It tests each contiguous candidate merge: for a given positive count the attainable mean is precisely the interval of its panel means, and attainable count-excess pairs are a convex cone. Consequently feasibility at the two endpoints implies feasibility for their mixtures. The two accepted q meshes agree on a 13-parent partition (four removed parent boundaries, all kept response bins are unchanged). At S4 reference weights 0, .25, .50, .75 and 1, the resulting joint absolute optima are respectively 1.46341e-4, 1.46193e-4, 1.46123e-4, 1.46155e-4 and 1.46187e-4. Paired-mesh defect 2.1e-13. The small nonzero finite joint residual survives this choice of common reference; the large original own-moment separation is not recovered by changing the shared reference.

The fallback chooses a maximal-count grouping, not an optimized physical scalar grid; it deliberately discards four independent moment constraints. The result neither replaces the evolved state nor estimates what an uncomputed finer panel would return. No directed rounding.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar; d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_cold_shared_reference_scan.py
python -W error early_pole_endpoint_cold_shared_reference_scan.py \
 --left "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --right "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{2,4}-x05-state-$d.json.gz \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
```
