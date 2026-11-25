# Earlier production onset: abundance-matched hot-tail reduction — 2026-11-25

Use the same energy-closed background model, g_bath=106.75 and chemical relaxation/H=100, initially empty daughters, equilibrium parent. Vary x0=M/T_initial and retune physical width Gamma/H(M) to approximate omega_phi h²=.12 without renormalizing any spectrum.

| x0 | width parameter | omega h² | p/T>.1 number fraction | max f (grid diagnostic) |
|---:|---:|---:|---:|---:|
| .1 (prior 480) | .1383 | .12015 | 43.59% | ~2e8 |
| .05 | .095218 | .11995 | 30.84% | 3.48e10 |
| .02 | .0553 | .12077 | 18.69% | 1.88e13 |
| .01 | .0353 | .12062 | 12.34% | 1.77e15 |
| .005 (480) | .0217 | .11997 | 7.85% | 1.51e17 |
| .002 (320 exploratory) | .0109 | .11922 | 4.04% | 4.74e19 |

For x0=.005, increasing from 320 to 480 cells while lowering qmin by a factor of 10 changes abundance by 0.089%, mean momentum by 0.144% and hot fraction by 0.0114 percentage points. Normalized continuity residual at 480 is 5.43e-9. Thus early onset is a much stronger lever than finite vacuum scalar mass at electroweak parent scales. This does not prove a realizable high-temperature scenario; the enormous peak occupancy, thermal masses, reheating history and interaction rates are unresolved. The x0=.002 case requires independent IR convergence. Local CLASS export quadrature/moment checks succeed even for the extremely narrow spectra, but perturbations are not yet computed for them.

[Reproduction script](reheating_matched_cases.py), [selected actual local outputs](results/reheating-matched-2026-11-25.json). Next compute synchronous CLASS spectra on x0=.005 independent grids, verify CDM hybrid and long-wave normalization, and explore whether a physical thermal mass cuts off the cold birth channel. Do not interpret threshold p/T>.1 as observationally constrained hot DM.
