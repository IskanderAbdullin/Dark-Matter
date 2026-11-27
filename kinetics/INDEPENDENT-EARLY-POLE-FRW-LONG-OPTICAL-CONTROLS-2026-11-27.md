# Two long controls for the conditional full H clock

Compare the [96-step expanding full-H run](INDEPENDENT-EARLY-POLE-FRW-LONG-EXPANDING-CLOCKS-2026-11-27.md) with two independent calculations using the initial dark-to-visible ratio .05135770288767686. Keep the same fixed dielectric responses, domain, native occupations and macro count. The first control follows exactly the same redshift path, multiplying all collision rates by .9752697401968827. The second has normal radiation-only collision rates, but ends at the shorter log-scale interval 7.021942129417556e−5.

| Model | Terminal scalar gain | Gain difference vs full H | Maximum count history difference vs full H |
| --- | ---: | ---: | ---: |
| Full H, identical path reference | 7.631411815e−6 | 0 | 0 |
| Same path with initial constant multiplier | 7.631400988e−6 | −1.083e−11 | 1.083e−11 |
| Shorter radiation-clock expansion | 7.631410208e−6 | −1.607e−12 | 1.367e−11 |

Compare corresponding macro indices; the shortened run does not sample identical redshift. Its smaller terminal discrepancy is not evidence that it is the more faithful dynamic approximation. The maximum count differences in both controls are under two parts per million of the full-H final gain, whereas the unrescaled radiation clock changes it by 2.668%. The final conditional rate multiplier evaluated from the energy ledger is .9752697398324044, differing from its initial value by 3.645e−10. Main and control trajectories use different midpoint bisection locations (full 80/83/91, same-path 91/96, shorter 37); this experiment does not resolve a separate dynamic-clock correction below that numerical uncertainty. The energy and entropy closure gates remain unchanged, and no late-time extrapolation follows.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
python -W error early_pole_joined_screened_frw_reference.py fine \
 --duration 7.2e-5 --steps 96 --expands --clock-scale .9752697401968827
python -W error early_pole_joined_screened_frw_reference.py fine \
 --duration 7.021942129417556e-5 --steps 96 --expands
```
A separate checkpoint path for each run is optional. Raw and comparative diagnostics: `results/early-pole-joined-frw-long-optical-controls-2026-11-27.json.gz`.
