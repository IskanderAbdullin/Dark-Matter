# Endpoint evidence and limits at the current checkpoint

This independent endpoint fixture uses an evolved scalar histogram, a chosen optical pole and accepted quadratures for restricted reaction bands. It does not determine the all-momentum, all-channel production history. The thermal CLASS and toy cooler-distribution benchmarks on the base branch are separate calculations.

| Question | Accepted evidence | What it leaves open |
| --- | --- | --- |
| Does increasing the scalar spatial resolution reduce the held cold response shape difference? | [Complete S16 two-mesh check](EARLY_POLE_ENDPOINT_COLD_S16_SPATIAL_PROGRESS.md). The three successive conditional forward TV distances on the same 17 bins are 2.023%, 0.451%, and 0.111%. | No asymptotic rate, no high-q tail bound, no extrapolation beyond the retained cold band. |
| Is the change carried by the explicit scalar factor or the bare kernel? | [One exact spatial replacement order](EARLY_POLE_ENDPOINT_COLD_SPATIAL_POPULATION_PATH.md), with the opposite order explicitly unresolved where native panels cross scalar edges. The [same-chart final-step clock control](EARLY_POLE_ENDPOINT_COLD_CLOCK_BAND_FACTORIZATION.md) isolates a different axis. | No order-independent spatial attribution until both crossed states use measured child coefficients. The bare kernel still combines multiple physical factors. |
| Is the warm pressure robust to quadrature order and the final clock? | Actual [G5 full history](EARLY_POLE_ENDPOINT_WARM_GAUSS5_HISTORY.md), [G5 final-clock replay](EARLY_POLE_ENDPOINT_WARM_GAUSS5_PRESSURE_CLOCK.md) and [G6 paired startup](EARLY_POLE_ENDPOINT_WARM_GAUSS6_STARTUP.md). Warm pressure may move even when total scalar number barely changes. | Full G6 history and cumulative temporal sensitivity; fixed-count and fixed-energy pressure intervals are conditional on already evolved moments. |

## Local verification in one standalone tree

At the PR revision `c818e2bee41a9b6bf7a8cb24f078003057fda616`, verification from a fresh unpacked checkpoint used only the published Python files and result inputs, with no other work directory on the import path. All 74 shipped endpoint unit tests and syntax compilation passed. Recomputed spatial scalar-factor bounds, S8/S16 native shared-reference responses, and G6 startup contrast each byte-compared against their committed receipts. The individual reproduction commands are in the notes linked above.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=. \
  python -W error -m unittest discover -p 'test_early_pole_endpoint_*.py'
python -W error -m compileall -q .
```

This is a local cloud test, not approval of the GitHub workflow. It does not validate a relic abundance or an Einstein–Boltzmann input spectrum. Before physical transfer inference the candidate needs a physical clock and momentum scale, full production/cooling and inverse-channel coverage with matched background with conserved energy and number budgets, an end-of-production extrapolation and high-q bound, and checks on the normalization and perturbation handoff. Agreement of finite held profiles, or of a normalized ballistic toy, cannot replace those steps.
