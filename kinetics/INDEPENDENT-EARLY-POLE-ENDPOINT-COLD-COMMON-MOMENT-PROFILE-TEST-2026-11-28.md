# A counterfactual with common scalar moments nearly matches the two cold profiles

When each chart keeps **its own** original fine-cell amounts and means, the closest S2 and S4 normalized forward-number profiles on the common q partition differ by about 2.02256% TV. This separate test changes the constraints: use the S2 coarser scalar cells and their S2 in-band counts and mean excess as the **same prescribed observations for both**. Use the S2 reference mass in both mean-excess constraints. Release the original S4 finer-cell moments, but keep each chart's own original charge and optical response. Each finite q panel then has a nonnegative scalar amount, subject only to these shared coarser constraints.

The jointly optimized normalized profiles become indistinguishable at the numerical precision of this finite test. With counts alone the nonnegative solver objective is about (1.4	imes10^{-17}) (direct witness (4.3	imes10^{-14})); with count and mean excess it is (3.68	imes10^{-12}), with an independent direct witness at the same level. Both valid q meshes give essentially the same results. These are finite witness residuals, not mathematical proof of exact equality.

The test shows that a large normalized separation is **not forced by the two grids and their fixed optical kernels alone** if one replaces the different retained scalar-moment constraints with common observations. The interval and mesh checks for the actual independent moments remain valid. This intentionally releases S4 information, fits a counterfactual profile and does not alter or explain the original evolution or absolute rates. In particular it is not justification to replace an evolved S4 state with the optimized profile or a proof of physical equivalence.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar; d=2026-11-28
python -W error early_pole_endpoint_cold_common_moment_profile_test.py \
 --left "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --right "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{2,4}-x05-state-$d.json.gz \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
```
