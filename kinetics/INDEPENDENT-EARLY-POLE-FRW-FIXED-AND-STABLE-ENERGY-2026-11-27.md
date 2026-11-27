# Fixed-scale control and a better-conditioned energy increment

The fixed-scale fine, 96-step control at exposure 7.2e−5 has gain 7.839594925e−6. Compared with the [expanding radiation reference](INDEPENDENT-EARLY-POLE-FRW-LONG-EXPANDING-CLOCKS-2026-11-27.md) at the same macro grid, prescribed redshift alone increases the gain by 0.012775%. This small difference is separate from the 2.668% change from including modeled dark energy in the expanding H clock. The fixed/radiation terminal interval mean rates are .1176803/.1177059. No late abundance or freeze-out is inferred.

Repeat the fixed-scale remainder after macro step 66 using an algebraically equivalent finite-energy difference. The reference bare energy is affine in the retained two moments, so its finite increment is evaluated directly in those moments. The smaller non-affine redshift increment uses a rationalized square-root difference and the FD occupation increment uses expm1 for small logit changes. The reciprocal gauge interaction is still included. Independent finite total-energy and redshift-trace fixtures agree in both archives.

The original fixed reference required bisections at macro steps 70 and 84. The switched stable-energy control reaches the end without bisection, has energy and entropy ledger defects +4.71e−13 and −2.79e−12, and changes the integrated gain by −1.424e−11 relative to the original fixed control (about 1.4% of the small prescribed redshift gain difference). A stage-matched single macro replay reproduces the stable implementation's step 67 diagnostics exactly. These results do not establish that the alternate energy algebra prevents every nonlinear convergence issue or that the redshift difference is temporally converged.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
python -W error early_pole_joined_screened_frw_reference.py fine --duration 7.2e-5 --steps 96
python -W error early_pole_joined_screened_frw_reference.py fine --duration 7.2e-5 --steps 96 \
  --stable-energy-from-step 67
python -W error early_pole_joined_screened_frw_reference.py fine --fixture
```
Raw: `results/early-pole-joined-frw-fixed-and-stable-energy-2026-11-27.json.gz`.
