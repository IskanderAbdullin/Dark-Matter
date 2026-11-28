# The completed S8 band reduces, but does not remove, cold shape differences

At the same D2 charge geometry and retained k=0.0015, integrate the actual independently evolved S8 endpoint x=.05 on both full q meshes. Each contains 100 panels on [1e-6,1.5e-5]; the nine ledger coefficients agree to 5.63e-14 relative between the two meshes. The fine S8 integrated forward-number response is 6.443407478e-7. Relative to S4 the forward-number, inverse-number, and forward-kinetic-excess integrals change by −0.08558%, −0.10000%, and −0.55514%. The full scalar yield changes by only −6.23e-7 in relative fraction. The two positive band-conditioned neutrality roots are about 2.913154e-9, not predictions of an ambient photon population.

For both spatial transitions use the **same exact 17-bin partition** common to S2, S4 and S8. Allow positive redistribution on each original q panel while retaining each chart's own scalar cell count and mean excess energy. Minimum achievable normalized TV:

| Channel | S2→S4 | S4→S8 |
| --- | ---: | ---: |
| Forward number | 2.02256% | 0.45149% |
| Inverse number | 2.02039% | 0.45099% |
| Forward kinetic excess | 1.87111% | 0.43506% |

The same q bin [3.30600e-6,3.97454e-6] remains the largest signed contributor in all six original comparisons. Assessing S4/S8 on their own finer 31-bin partition instead gives forward conditional TV 0.53775%, inverse 0.53296%, and forward excess 0.49429%. Thus the shape discrepancy decreases even on a finer display of that pair. These are comparisons between finite computed spectra at held fitted poles, not an asymptotic order, an S16 result, or a physical convergence claim.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
base=results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
python -W error early_pole_endpoint_cold_common_grid_spatial_ladder.py \
 --s2 "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --s4 "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --s8 "$p"8-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{2,4,8}-x05-state-$d.json.gz --archive "$base"
for obs in forward_number inverse_number forward_excess; do
 python -W error early_pole_endpoint_cold_profile_exact_finite_distance.py \
  --left "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz --right "$p"8-endpoint-continuum-cold-{1,2}.jsonl.gz \
  --states "$s"{4,8}-x05-state-$d.json.gz --archive "$base" --levels 4 8 --observable "$obs"
done
```
