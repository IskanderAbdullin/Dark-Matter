# Common cold scalar moments almost, but not exactly, align absolute responses

The earlier common-S2-moment control allowed every forward profile to be normalized separately. This test instead keeps the finite absolute response in each exactly shared q bin. It uses one nonnegative panel-occupation reweight on each chart to predict forward number, inverse number, and forward kinetic excess together. The positively weighted loss sums three L1 discrepancies; each field is divided by its own *original S2 band total*, not by a separately fitted total on each chart. No gain or amplitude is freely fit. Forward spontaneous terms remain fixed; inverse terms have no offset.

| Prescribed in-band scalar moments | Forward number alone | Inverse number alone | Forward excess alone | Joint loss |
| --- | ---: | ---: | ---: | ---: |
| Own chart, count | 0.01587750 | 0.01285430 | 0.00380504 | 0.07255557 |
| Own chart, count + mean excess | 0.04012591 | 0.04002355 | 0.04619394 | 0.12645695 |
| Shared S2, count | numerical zero | numerical zero | numerical zero | 0.00002377917 |
| Shared S2, count + mean excess | 0.00003320908 | 0.00002832751 | 0.00000105986 | 0.00014757540 |

In the shared count+mean variant, the three discrepancies for the **same witness** are 0.00007725791, 0.00006925762 and 0.00000105986. Thus the previous near-zero *normalized* forward-profile witness does not imply joint amplitude equality. Still, nearly all the original absolute separation can be removed in this counterfactual. That conclusion requires releasing the original finer S4 scalar moments. When the charts keep their own recorded moments, their joint optimum remains near the original loss (0.12757148); separate per-field estimates cannot be combined as a common witness. At least for the count-only case this distinction matters noticeably, including in the same-moment counterfactual.

Chart optical fits, charge response and endpoint remained fixed. Both accepted q meshes reproduce the optima to 1.9e-13 absolute. Explicit primal, dual and direct witness residuals are checked; no directed rounding, continuum inference, physical abundance calculation, or claim about scalar phase space beyond the archived band.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar; d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_cold_joint_absolute_profile.py
python -W error early_pole_endpoint_cold_joint_absolute_profile.py \
 --left "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --right "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{2,4}-x05-state-$d.json.gz \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
```

## Independent finite loss bounds

The signed support weights are clipped to [-1, 1]. For any such weights the signed projection is no larger than the L1 norm; minimizing it separately over every stored scalar parent gives a computed lower bound without relying on approximate LP dual feasibility. For the upper witness the solver weights are clipped nonnegative and explicitly projected onto their separate simplices. In the shared count+mean control the two bounds are 0.0001475753957192 and 0.0001475753957561, with no directed rounding. The own count+mean bounds are 0.1264569558607 and 0.1264569560867. The finite interval [3.306e-6, 3.97454e-6] accounts for 57.06% of the native and 46.92% of the counterfactual joint residual *at these projected optimum witnesses*. The split need not be unique across different optimum witnesses. Exact signed weight and absolute residual vectors are in the JSON.
