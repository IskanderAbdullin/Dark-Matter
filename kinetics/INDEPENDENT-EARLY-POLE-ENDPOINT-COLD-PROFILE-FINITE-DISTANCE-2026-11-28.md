# Minimum normalized cold separation over both finite occupation families

The separate-share lower bound from the [common profile envelope](INDEPENDENT-EARLY-POLE-ENDPOINT-COLD-NORMALIZED-PROFILE-ENVELOPE-2026-11-28.md) can understate the actual minimum: a chart may not attain several individual bin extrema at once. In the same S2/S4 finite deformations, now optimize the two **entire normalized profiles simultaneously**. Keep the exact common q groups and independently fixed charge and pole for each chart.

For each original scalar cell, let its chosen nonnegative occupation be a convex mixture of its count or count-and-mean vertices. Rescale all vertex weights for a given chart by one inverse total forward coefficient. After this change of variables, the original shell-weight sums equal the chart's rescaling variable, and the rescaled total equals one. All common-bin weights are now linear in the variables. A third variable per common group bounds its absolute difference, giving a linear program to minimize half their sum.

| Moment constraints | Previous separate-share TV lower bound | Minimum joint finite TV | Original step-shape TV |
| --- | ---: | ---: | ---: |
| Count in each cell | 0.729141% | **0.883307%** | 2.024020% |
| Count and mean excess | 2.022550% | **2.022556%** | 2.024020% |

The optimized chart profiles witness those minima on the declared variables. Separate direct variation, primal feasibility, and solver dual checks close within (3.1	imes10^{-10}) absolute fraction. The optimum differences between the accepted panel meshes are below (3.6	imes10^{-14}). These are numerically resolved finite optima, not interval-arithmetic certificates, pressure or cosmological errors, nor a solution of the continuous in-cell reweighting problem.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_cold_profile_exact_finite_distance
p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar; d=2026-11-28
python -W error early_pole_endpoint_cold_profile_exact_finite_distance.py \
 --left "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --right "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{2,4}-x05-state-$d.json.gz --levels 2 4 \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
```

## Independent inverse and excess projections

Repeat the very same normalized finite minimization for inverse number and forward excess. The inverse has no spontaneous offset, so verify a strictly positive lower total from the feasible nonnegative per-shell vertices before normalizing. The forward excess retains its spontaneous term.

| Normalized observable | Minimum with counts only | Minimum with counts and means | Original step-shape TV |
| --- | ---: | ---: | ---: |
| Forward number | 0.883307% | 2.022556% | 2.024020% |
| Inverse number | 0.796710% | 2.020395% | 2.022563% |
| Forward excess | 0.151856% | 1.871110% | 1.907637% |

Thus count-only gives less protection for the **shape** of forward excess than for number. A known in-band mean substantially constrains all three original profile differences, but the norms are distinct diagnostics, not physical errors. The forward-number receipt is unchanged. Use the preceding CLI with `--observable inverse_number` or `--observable forward_excess` to reproduce the other two.
