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

## Checks that can change the interpretation

1. Compare the local shift at matched held-screening split widths
   against the simultaneous trajectory, including its *history* and
   not just the terminal scalar total.
2. Compare time widths 16 and 32 at both the original hard aperture
   and the fourfold aperture. Only a crossed comparison identifies
   interaction between support and time width in this finite model.
3. Compare radial orders one, two and four with the hard support fixed.
   Keep that narrower diagnostic distinct from transfer/angle
   quadrature or an extrapolated cold spectrum.
4. Test all four daughter preparations on one simultaneous graph to
   check nonadditivity. Test the local perturbation while matching
   total retained daughter number and bare energy. The [static
   two-moment control](INDEPENDENT-EARLY-POLE-DAUGHTER-COMPLEMENT-RETUNE-2026-11-27.md)
   already shows that the initial local suppression is unchanged to about
   33 ppm when the complement repairs the number and bare energy.
   Its finite-time result remains a separate check; the two matched
   moments do not fix the dressed energy or the full phase-space shape.

The [longer warm-marker redshift](INDEPENDENT-EARLY-POLE-LIVE-WARM-REDSHIFT-STRICT-2026-11-27.md)
is a separate prescribed expansion experiment. Neither it nor the
short held-response results establishes the abundance, perturbation
transfer function, observational viability or halo properties.
