# Early-pole first-window controls: what survives refinement

This note links numerical controls from the same prepared fine hard
state and short exposure (3.3333333333333335e-7), rather than inferring
a physical production rate. The native preparation and the continuous
cell-conservative preparation are different archived inputs. A
swap makes one controlled hybrid of those inputs on common support;
it is not a different species, interaction or microscopic model.

## Stable and unstable conclusions

- **A simultaneous time integrator yields a resolved finite contrast.**
  In the held dielectric model, the two preparations differ by about
  1.467% in mean cold gain relative to the native result. A tenfold
  tolerance tightening changes the two means by less than 0.5 ppb
  within that one chart and discretization. It does not tighten the
  error due to hard support, elastic packet order or the physical model.
  See the [unsplit comparison](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-UNSPLIT-CHART-2026-11-27.md)
  and [tolerance test](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-UNSPLIT-TOLERANCES-2026-11-27.md).
- **The local daughter contribution changes sign between time solvers.**
  Its mean shift is about -1.00e-5 with 16 split steps and +5.17e-4
  with simultaneous integration. Holding the response does not remove
  the difference at 16 steps. At that width, a small scalar total
  conceals a 91.7% cancellation among opposing cold-bin effects.
  See the [held-response test](INDEPENDENT-EARLY-POLE-DAUGHTER-HELD-RESPONSE-CONTROL-2026-11-27.md)
  and [mode audit](INDEPENDENT-EARLY-POLE-DAUGHTER-COLD-MODE-EFFECT-2026-11-27.md).
- **The outside daughter reduction is more robust across the tested
  solvers.** It is about -7.42e-4 in the two-channel split control
  and -7.72e-4 when evolved simultaneously with fixed dielectric
  weights. Unlike-only scattering recovers the two-channel split
  effect within 1.42%; the two elastic contributions are not additive.
  This is not a universal partition of a physical production rate.
  See the [elastic control](INDEPENDENT-EARLY-POLE-OUTER-DAUGHTER-ELASTIC-ROLES-2026-11-27.md)
  and [prepared-input unsplit results](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-PREPARED-UNSPLIT-2026-11-27.md).
- **The hard support is not converged.** One, two and four subdivisions
  of the original active hard cells give native/continuous contrasts
  of 2.470%, 2.330% and 2.578% at the same split width. The reversal
  occurs even though the corresponding static continuous source
  changes much less after depth two. Constant closure gates are not
  a test of momentum-resolution convergence.
  See [hard aperture depth](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-APERTURE-DEPTH4-2026-11-27.md).
- **Elastic quadrature moves the finite comparison at fixed support.**
  Both radial packet orders doubled changes the same split contrast
  from 2.470% to 2.335%. This is a separate control from hard support,
  not an alternate estimate of its error.
  See [elastic radial quadrature](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-FINITE-ELASTIC-REFINEMENT-2026-11-27.md).

## Results that narrow the numerical and preparation uncertainties

| Question | Matched control and result |
|---|---|
| Is the local sign different with a smaller time step? | [64 held split steps](INDEPENDENT-EARLY-POLE-LOCAL-DAUGHTER-HELD-TIME64-2026-11-27.md) yield +3.85e-4 rather than the 16-step -1.00e-5. The simultaneous result is +5.17e-4, and the [64-step mode audit](INDEPENDENT-EARLY-POLE-DAUGHTER-COLD-MODE-EFFECT-2026-11-27.md) finds all eight terminal differences positive. |
| What changes when only the elastic substep is refined? | [16 macro × 4 elastic](INDEPENDENT-EARLY-POLE-LOCAL-DAUGHTER-ELASTIC-SUBCYCLE-2026-11-27.md) gives +2.13e-4; 64 × 1 gives +3.85e-4 at the exact same elastic width. Both the inner step and outer interleaving matter. |
| What does the intermediate outer width show? | [32 macro × 2 elastic](INDEPENDENT-EARLY-POLE-LOCAL-DAUGHTER-FIXED-INNER-MIDDLE-2026-11-27.md) gives +3.47e-4 with the same inner width as 16 × 4 and 64 × 1. The observed difference ratio is 4.53 rather than the single outer-quadratic value 5. Separate quadratic diagnostics leave about 1.20e-4 relative to the simultaneous result. |
| Can a time correction measured on one hard mesh be reused at another? | The [16/32 × depth 1/4 split test](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-TIME-APERTURE-CROSS-2026-11-27.md) finds the continuous preparation's time-refinement gain is about 49% smaller at depth four. |
| Does the latest radial refinement remain large? | Increasing both packet radial orders from 2 to 4 moves the two rates by -168 and -129 ppm at the same steps and support, much less than the 1-to-2 movement. [Radial depth four](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-ELASTIC-RADIAL-DEPTH4-2026-11-27.md) is not a transfer, angle or support convergence test. |
| Do local and outside daughter changes add? | [Full held-chart substitution](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-UNSPLIT-DAUGHTER-ADDITIVITY-2026-11-27.md) finds the positive local effect offsets about 67% of the negative outside effect; the small residual is [unchanged by tolerance tightening](INDEPENDENT-EARLY-POLE-PREPARED-CHART-TOLERANCE-2026-11-27.md). This is not the elastic-channel interaction. |
| Is the local rebound robust to matching retained bare moments? | [Two-moment compensation](INDEPENDENT-EARLY-POLE-LOCAL-DAUGHTER-RETUNED-TIME-2026-11-27.md) changes the initial suppression only 33 ppm but reduces the positive finite mean shift by about 83%. The [same outside compensation alone](INDEPENDENT-EARLY-POLE-DAUGHTER-COMPLEMENT-FACTORIAL-2026-11-27.md) explains almost all that reduction for this particular retuning. |

The outstanding checks include reducing the inner step at 64 macro steps,
the simultaneous chart on the fourfold aperture, and independent
compensation controls. Neither a small integrator
tolerance sensitivity nor exact matching of two bare moments
bounds the uncertainty from the support or physical kernel.

The [longer warm-marker redshift](INDEPENDENT-EARLY-POLE-LIVE-WARM-REDSHIFT-STRICT-2026-11-27.md)
is a separate prescribed expansion experiment. Neither it nor the
short held-response results establishes the abundance, perturbation
transfer function, observational viability or halo properties.
