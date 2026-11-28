# Separate the local energy threshold from FRW time-step refinement

Use the same two prescribed cosmologies, the same two warm-marker states
in each, and the same initialized comoving cold and hard sectors. The
72-step profiles differ only in the elastic energy threshold. The
72- and 144-step strict profiles share the 5e-15 threshold. The receipt
cross-checks initial count rates and the initialized gap, audits the
redshift trace subtraction, and applies the original global gates to
all four arms of each profile.

| Profile | Warm expanding-minus-fixed (ppm) | Cold expanding-minus-fixed (ppm) | Largest trace-subtracted energy defect |
|---|---:|---:|---:|
| Original 72 | 36.020441496 | 35.917497132 | 1.8261e-12 |
| Strict 72 | 36.020441998 | 35.917498009 | 1.7764e-15 |
| Strict 144 | 36.020607174 | 35.899086954 | 4.4538e-15 |

At 72 macro steps the tighter gate reduces the largest residual
roughly one thousandfold while moving all reported absolute gains
by less than 8e-6 ppm. In contrast, doubling the steps at the same
strict threshold changes the absolute warm gains by about 7.68 ppm
and cold gains by about 30.92 ppm. The paired warm redshift contrast
changes only about 0.000165 ppm. The benefit from the threshold is
in the ledger, while the tested absolute-gain change mostly comes
from the time grid at the tighter threshold.

The first comparison isolates the numerical energy threshold at fixed
macro width, while the second isolates the macro-width change at fixed
threshold. These are three corners, not an estimate of interaction.
The warm marker is on a prescribed radiation clock and held-polarization
model. The small relative redshift pair is not a physical relic
prediction and the data do not establish formal convergence.

Use the native fixed/expanding raw pairs from the original 72-step,
strict 72-step and strict 144-step runs with
`early_pole_warm_frw_isobudget_comparison.py`. Each argument takes
fixed, then expanding raw JSON.
