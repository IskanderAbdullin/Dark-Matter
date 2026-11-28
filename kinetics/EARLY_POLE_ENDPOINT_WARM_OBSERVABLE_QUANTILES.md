# Where the warm-grid effect lives in momentum

Compare the evolved fine S4 native-warm and W2 scalar histories (same charge chart d=2 and same cold scalar cells). For each of the four published snapshots, integrate each histogram, rather than taking cell-center samples. Hold the numerical mass at the native fit for the shape comparison, then reevaluate W2 at its own fit. The integration returns independent cumulative distributions for count, excess over rest energy, the numerator of pressure, and the scalar-count-weighted squared velocity. The pressure result intentionally omits common factors; the last observable is not the pressure. Use the packed endpoint states in `results/` to replay the diagnostics.

At x=.05:
- Scalar count changes by only -3.095e-7 (relative). Its median q is about 3.55708e-6 for either chart. Above q=.01, native and W2 each have fewer than ten particles per million, but over 99.8% of the pressure numerator.
- The median pressure-weighted q shifts from 1.84969 to 1.83469. The pressure numerator changes by -6.23477%; reevaluating W2 at its own scalar mass changes it by a further -1.33e-11 relative. The large relative response is not explained by the fitted mass.
- The velocity-squared diagnostic places the median at 4.2283e-6 and 4.2082e-6, respectively. It changes by -1.2472%, not -6.2348%. The high tail carries 34.38% and 33.55% of this diagnostic; velocity-squared saturates at large q.
- Over x=.01, .02, .04 and .05, pressure changes between -6.32% and -6.23%, while this velocity-squared diagnostic changes between -1.25% and -.14%. These are frozen snapshots of paired evolved trajectories, not a continuous-in-time bound.

From the `kinetics/` project directory with the endpoint harness dependencies available:
```sh
python -W error -m unittest test_early_pole_endpoint_warm_observable_quantiles
for i in 01 02 04 05; do
  python -W error early_pole_endpoint_warm_observable_quantiles.py \
    results/early-pole-fine-charge2-scalar4-x${i}-state-2026-11-28.json.gz \
    results/early-pole-fine-charge2-scalar4-warm2-x${i}-state-2026-11-28.json.gz --x 0.${i}
done
```

The count integral is checked analytically, all integrals and reported quantiles use a 24/48-node comparison, and each quantile is checked by integrating its remaining upper tail separately. A massless fixture supplies exact moments and quantiles; it does not alter the positive fitted masses in the receipts. This is sensitivity within two numerical reconstructions of the same restricted collision model. It is not a CLASS input or a validated physical transfer prediction and does not establish either chart as converged.
