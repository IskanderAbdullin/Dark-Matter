# One-node intervention on the held warm source

At x=.01 the largest directly measured G8→G9 change selects the warm original cell [3.96368,5.72884]. Upgrade that cell only, retaining G8 elsewhere (300 total nodes rather than 299 for G8 and 321 for G9). As a matched-budget control, use the next large parent [5.72884,8.28009] instead. Test both frozen choices unchanged at x=.05. The newly assembled direct source reproduces the expected single-parent substitution, including independent forward and inverse ledgers.

The trained parent captures 85.03% of the signed whole G8→G9 change at x=.01 and 68.17% at the held x=.05, versus 3.05% and 2.15% for the equally expensive control. At x=.05 the trained aggregate residual against Volume64 becomes -0.0034105%, versus -0.0275876% with G8. Its sum of absolute cell residuals becomes 0.0206304%, versus 0.0286017% with G8 and 0.0178374% with G9. The very small trained aggregate difference at x=.01 does not bound its cellwise error.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results
python -W error early_pole_endpoint_warm_single_g9_upgrade.py \
 --states "$p/early-pole-fine-charge2-scalar4-x"{01,05}"-state-$d.json.gz" \
 --refs "$p/early-pole-endpoint-warm-source-cells-x"{01,05}-$d.json" \
 --lock "$p/early-pole-warm-gauss9-10-12-frozen-source-$d.json"
```

Both clocks are held states from the same baseline history. This tests a quadrature intervention, not whether its benefit survives evolution or holds for arbitrary inputs.
