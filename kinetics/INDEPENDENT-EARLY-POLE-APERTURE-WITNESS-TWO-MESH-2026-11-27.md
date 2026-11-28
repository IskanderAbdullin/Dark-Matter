# The dressed sign for the strongest witness depends on resolving scalar momentum

Continue the [frozen reconstructed aperture](INDEPENDENT-EARLY-POLE-APERTURE-WITNESS-REPROJECTION-2026-11-27.md) with 16 and 32 scalar bins and 16 and 32 hard subdivisions independently. The same original parent-only bounded witness is used, with matching archived native moments before projection. No distributions are evolved.

| Grid | Witness fraction | Dressed rate, 8 scalar / 16 hard | 16 scalar / 32 hard | 32 scalar / 32 hard |
| --- | ---: | ---: | ---: | ---: |
| Fine | 85% | +.010028 | +.000354 | +.000787 |
| Fine | 97% | +.000328 | −.010823 | −.010308 |
| Third | 85% | +.035483 | +.031130 | +.029589 |
| Third | 97% | +.030451 | +.025386 | +.023515 |

The fine 97% case becomes negative even under the dressed aperture at both finer resolutions and at both hard subdivisions. Fine 85% remains slightly positive in the four dressed refinements, and third stays positive. This remains an initial closure-model comparison, not a trajectory, and the positive fine 85% residual is sensitive relative to its own magnitude. No inference about full dynamical screening, finite widths, or freeze-out follows.

A separate bare integration uses the original step density in scalar phase volume per native shell and integrates in (q^3), splitting at all native density discontinuities and at a 64-bin auxiliary mesh. Sixteen versus 32 Gauss nodes per panel differ by less than 7e−7, and 512 geometric volume-center bins agree within 3.2e−6 absolute. Bare 85% rates from the finer independent quadrature are −.00465563 on fine and −.00767136 on third; both are already negative at 16 volume-center bins. The apparent fine bare reversal loss at eight bins therefore comes from the coarse mode sampling in this model, not the conserved scalar number under rebinning.

From `kinetics/`:

```sh
OPENBLAS_NUM_THREADS=1 python early_pole_scalar_grid_aperture_witness_reference.py
for g in fine third; do
  for n in 16 32; do
    OPENBLAS_NUM_THREADS=1 python early_pole_live_aperture_witness_reprojection_reference.py \
      "$g" --subdivisions 16 32 --scalar-bins "$n"
  done
done
```

Raw: `results/early-pole-aperture-witness-two-mesh-2026-11-27.json.gz`.
