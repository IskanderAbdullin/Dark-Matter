# A triangle-compatible live scalar shape sensitivity test

Apply the [triangle eligibility projection](INDEPENDENT-EARLY-POLE-SCALAR-TRIANGLE-SECANT-2026-11-27.md) to the continuous-time, isolated scalar shape evolution, keeping all archived charged inputs fixed. For the restricted graph, accepted outgoing pairs satisfy the incoming-momentum triangle and their interpolation mean conserves the instantaneous Hartree energy gradient. Use a positive mobility from the ratio of the usual held-stencil event flux to its affinity, then use the projected entropy gradient. This is a numerical proxy and not a recalculation of the screened matrix element at new scattering angles. The regular graph serves as a matched control with the same fitted initial moments.

| Archive and outgoing interpolation | Final scalar shape L1 / initial count | Final native source | Final 97% witness source |
| --- | ---: | ---: | ---: |
| Fine, regular | .85877 | +.078646 | −.0046072 |
| Fine, triangle adjacent-only | .83135 | +.078668 | −.0045701 |
| Fine, triangle + one neighbor | .84210 | +.078659 | −.0045901 |
| Third, regular | .70547 | +.070042 | +.039343 |
| Third, triangle adjacent-only | .67708 | +.070193 | +.040904 |
| Third, triangle + one neighbor | .68612 | +.070138 | +.040319 |

No witness sign changes relative to the regular graph on these sixteen-cell, order-six controls at exposure 2e−7. The difference in final witness source reaches about 0.8% on fine and 4.0% on third for adjacent-only projection. This is sensitivity to the whole projected weak event construction, not simply a measure of the small number of rejected packets. Positive fixed-count coordinates give total scalar number defects below 3.5e−18 absolute and evaluated Hartree energy changes below 6.8e−21 across these numerical fixtures. No change to the physical coupled solver or freeze-out inference follows from these checks.

From `kinetics/`:

```sh
for g in fine third; do
  OPENBLAS_NUM_THREADS=1 python early_pole_scalar_triangle_live_proxy_reference.py "$g"
done
```

Raw: `results/early-pole-scalar-triangle-live-proxy-2026-11-27.json.gz`.
