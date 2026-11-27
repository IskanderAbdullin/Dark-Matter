# The prepared same-flavor relaxation tangent is dominated by small transfers in the tested sector

Starting from the fine and third moment-matched 97% witness or their untouched native steps, hold the screened same-flavor t-pole response fixed. The frozen collision increments conserve each flavor's number and close its variational energy ledger. Evaluate the scalar source at symmetric displacements along that collision vector. Repeat with both a finer Q/angle quadrature and three upper bounds on the plasma-frame momentum transfer, keeping the lower bound at 1e−4.

| Upper Q | Fine witness tangent | Third witness tangent |
| --- | ---: | ---: |
| .04 | +4.1129e7 | +8.2557e6 |
| .08 | +4.1722e7 | +8.2364e6 |
| .16 | +4.1885e7 | +8.2421e6 |

These are directional derivatives per squared reference exposure, not source values. The fine witness starts negative, while the third witness at this reconstruction starts positive. The untouched step tangent is negative on both archives for every tested cut. Reducing the cut from .08 to .04 changes the fine witness tangent by about 1.4%; widening it from .08 to .16 changes it by about 0.4%. The third changes are each below 0.25% on the finer quadrature. In the fine .08 test, the parent-only tangent is +4.1868e7 and the daughter-only tangent is −1.4580e5. The parent portion also dominates on third. These calculations independently support the parent aperture mechanism observed in the [matched short evolution](INDEPENDENT-EARLY-POLE-WITNESS-ELASTIC-SECTOR-ATTRIBUTION-2026-11-27.md).

The pole proxy is not a complete same-flavor collision operator; the finite higher-Q sample is not a matching prescription. No behavior below the artificial lower Q cut, other scattering channels, nonlinear tangent extrapolation to a crossing time, or populated physical witness is implied. Holding screening fixed is part of this fixture.

From `kinetics/` with `OPENBLAS_NUM_THREADS=1`:

```sh
for g in fine third; do
  for n in 0 1; do
    python early_pole_frozen_same_soft_cut_sensitivity_reference.py "$g" --level "$n"
  done
done
```

Raw: `results/early-pole-witness-frozen-same-soft-Q-cut-2026-11-27.json`.
