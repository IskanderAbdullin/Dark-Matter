# Quadratic velocity matching does not align the finite ballistic response

This is a separate stress test of the five **frozen** positive scalar distributions W1, W4 and G2–G4 at x=.05. It is motivated by the use of velocity and full phase-space distributions in [Ballesteros, Garcia and Pierre](https://arxiv.org/abs/2011.13458) and [Huang, Li and Yu](https://arxiv.org/abs/2306.00065). We do not use their calibration, reproduce their transfer functions, or map our dimensions to cosmological wavenumbers.

For an isotropic set of particles with frozen velocity magnitude v=q/sqrt(q²+M²), consider a density modulation in otherwise unforced free motion over a dimensionless phase b. Averaging exp(i b v cos theta) over orientations gives sinc(bv). Report the nonnegative loss D=1−mean[sinc(bv)], not a power spectrum. With common phase b the leading loss is b² times mean(v²)/6. In the separate matched-RMS control choose b separately per chart so that u=b sqrt(mean(v²)) is identical, making the quadratic coefficient exactly u²/6 on each chart. Its next coefficient is minus the standardized fourth velocity moment divided by 120.

At common W1 measurement mass M=0.0008731829968 and u=1:

| Comparison | G4 minus G3 loss | G4 minus W1 loss |
| --- | ---: | ---: |
| Same raw phase | −2.555e-9 | −4.062e-7 |
| Matched RMS phase | +1.100e-5 | +1.727e-3 |

In the RMS-matched case, the G4/G3 difference comes almost entirely from q<.01, because aligning the overall velocity variance changes how much phase accrues to this dominant cold population. In the same-raw-phase arm those phases are identical and the smaller warm-population change dominates the difference. Thus the larger matched-RMS effect is a feature of this calibration test, not a larger change at one fixed physical scale. At formal M=1 the u=1 matched G4/G3 difference instead remains small and negative, −2.408e-9; the sign is not universal. The massless limit reproduces the identical single-speed angular law for all charts.

Every accepted probe closes the numerical 96/192 within-cell integration check and the code checks the small-phase angular expansion independently. The heavier formal mass and largest phases include explicitly unresolved probes, which are excluded instead of plotted as accepted. At the highlighted M near .000873 and u=1, paired discrepancy is below 1.4e-17. No metric forcing, expansion, collision source, transfer function or background cosmological parameters enter this test.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_frozen_ballistic_shape
python -W error early_pole_endpoint_frozen_ballistic_shape.py \
 --native "$p"-x05-state-$d.json.gz --w4 "$p"-warm4-x05-state-$d.json.gz \
 --gauss "$p"-gauss{,3,4}warm-x05-state-$d.json.gz
```
