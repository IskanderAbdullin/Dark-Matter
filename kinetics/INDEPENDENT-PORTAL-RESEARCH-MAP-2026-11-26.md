# Independent portal calculations: what each result does and does not test

This is a navigation aid for the independent numerical work in
[PR #2](https://github.com/IskanderAbdullin/Dark-Matter/pull/2).
It does not replace the upstream research plan or the separate
energy-closed and BEST calculations in the repository. The
reproducible controls below are not observational constraints.

| Question | Evidence available | What is still missing |
| --- | --- | --- |
| Does the massive-shell projection conserve the represented moments? | Exact and independent checks of number, energy, pressure/trace, positivity, and unresolved cells; see the PR's projection and restart references. | Accuracy of a specific complete physical production model does not follow from those geometry checks. |
| Does this *low-reheat* spontaneous branch make a dominant 86 eV species? | The stable rest count is set near Ωh²=.120. The shifted first-order scalar rest count is only about 3.05×10⁻¹⁰ or 2.75×10⁻⁹ for the two prescribed histories. | The result is limited to these sources and rate assumptions. It cannot decide a branch using strong stimulated production, different portal parameters or new channels. |
| Is the scalar contrast robust to integrating the hierarchy? | A matched native fluid-closure test, the corresponding no-closure positive split, and a memory-bounded independent relativistic calculation on recorded metrics trace the discrepancy to the approximate grouping. | No global hierarchy convergence or coupled matter transfer from the independent fixed-metric run. |
| Can the *warm 50 MeV* reference be treated as numerically converged? | Restarted q16/q32 runs reach 10.442 MeV on matched grids. Their count difference is small, but absolute time refinement and higher spectral moments must be audited separately. | The short matched endpoint is not the 3.43 MeV long high-support run or a late relic. Earlier exterior phase distributions are not present in the restart archive. |
| Can purely ballistic forward depletion rescue the initial warm excess? | An isolated loss-only bound controls Bose gain even with all new scalar momenta and gives a stable lower proxy around 56 at the late anchor. | It does not control inverse repopulation, elastic redistribution, unknown inherited scalar tails or another physical history. |
| How sensitive are the hard thermal folds to the perturbative range? | Separate Born and strict-NLO extensions, exact measured-window subtraction, and an isolated toy high-mass pole test. | Full physical thermal QCD, unsupported spectral domains and threshold effects are not bounded. |

## Low-reheat branch: follow the count, then the spectrum

1. [The fixed-loss retarget](INDEPENDENT-FINITE-COOLING-FROZEN-COUNT-RETARGET-2026-11-26.md)
   shifts illustrative temperatures by 1.11–1.31 keV. The
   [actual re-evaluation](INDEPENDENT-FINITE-COOLING-SHIFTED-COUNT-2026-11-26.md)
   leaves only a 1.50–2.75 eV local tangent update. These are
   still conditional spontaneous forward counts, not a fitted
   or self-consistently iterated reheating temperature.
2. [Compare matched scalar spectra](INDEPENDENT-FINITE-COOLING-SHIFTED-SPECTRUM-2026-11-26.md)
   and [the independent relativistic hierarchy control](INDEPENDENT-FINITE-COOLING-RELATIVISTIC-LOS-2026-11-26.md).
   A scalar density contrast can differ strongly between histories
   while its physical count is much too small to dominate the dark
   rest density in either. A large relative scalar contrast alone
   is not a large matter effect.
3. [Combine both tracers and baryons on the same fixed metric](INDEPENDENT-FINITE-COOLING-STABLE-PAIRED-LOS-2026-11-26.md).
   The direct squared-contrast change is below 7×10⁻⁸ in absolute
   value at the five sampled modes, with no backreaction. The
   [separate homogeneous-clock bound](INDEPENDENT-FINITE-COOLING-SHIFTED-CLOCK-2026-11-26.md)
   concerns a converted-versus-unconverted proxy after 1 MeV,
   not perturbations, nucleosynthesis or measured expansion.

## Warm branch: preserve what the checkpoints actually contain

1. [The extended source support](INDEPENDENT-WARM-SIGNED-Q32-PRODUCTION-2026-11-26.md)
   preserves the inner nodes while adding outer ones. The
   [matched count and shape test](INDEPENDENT-WARM-EARLY-SUPPORT-TIME-INTERACTION-2026-11-26.md)
   shows why comparing mismatched time grids would mix numerical
   effects. The [sixteen-block follow-up](INDEPENDENT-WARM-EARLY-SUPPORT-16-BLOCKS-2026-11-26.md)
   verifies that the count support signs persist; higher moments
   remain a separate check. The [q16 time ladder](INDEPENDENT-WARM-EARLY-TIME-REFINEMENT-2026-11-26.md)
   is a fixed-support result.
2. [The atomic restart](INDEPENDENT-WARM-SERIAL-RESTART-2026-11-26.md)
   carries actual newly recorded exterior cohorts. For older cohorts
   known only by moments, [bound their omitted gravity](INDEPENDENT-WARM-EARLY-OMITTED-GRAVITY-2026-11-26.md)
   instead of inventing an occupation history; an instantaneous
   gravity bound does not bound integrated yields.
3. [The frozen split-order test](INDEPENDENT-WARM-FROZEN-SPLIT-ORDER-2026-11-26.md)
   finds that switching the implicit stepper alone does not remove
   the schedule sensitivity. This does not select a production
   integrator for the full expanding calculation.
4. [The closed-reservoir bound](INDEPENDENT-WARM-PAIR-ENERGY-AND-BAND-2026-11-26.md)
   allows redistribution but assumes a Bose cap. The stricter
   [parent-weighted gain](INDEPENDENT-WARM-PAIR-PARENT-WEIGHTED-BOSE-CEILING-2026-11-26.md)
   and [shell-resolved loss comparison](INDEPENDENT-WARM-PAIR-SHELL-SURVIVAL-2026-11-26.md)
   close such a cap only for ballistic forward-only evolution.
   Do not combine the two assumption sets to claim the cap is
   established in the general signed model.

The [hard thermal x=16 comparison](INDEPENDENT-HARD-VECTOR-NLO-X16-2026-11-26.md)
is a separate sensitivity test on specified mass bands and a
fixed radiation clock. Its conditional temperatures should not
be used as the low-reheat source temperatures above.

## A rate-independent number check

The [closed number budget](INDEPENDENT-PORTAL-CLOSED-NUMBER-BUDGET-2026-11-26.md)
applies across both archived branches **only** after specifying no
further change of total stable-plus-scalar number per visible entropy.
Even arbitrary equal-number conversion then cannot obtain a
dominant 86 eV rest component at the reference total .120.
This complements, rather than assumes, the narrower kinetic loss bound.
The [passive heavy-recycling exposure](INDEPENDENT-WARM-BALLISTIC-HEAVY-EXPOSURE-2026-11-26.md)
integrates positive formation from the warm 10.442 MeV anchor to
kinematic closure on those unchanged represented supports. It
credits only about 10⁻⁶ of the initial count, despite appreciable
relative sensitivity to the outer momentum cells. It does not
bound an evolving population with repopulation or older tails.

## Evidence boundary

Only the exact tested models, grids, and signed controls are
reported. Closing the physical question still requires a viable
production mechanism at the intended abundance, a self-consistent
transfer test of that mechanism, and only then an observational
comparison. Distinguish a failed numerical approximation, a
conditional exclusion of one mechanism, and lack of evidence
about the general model. The links above point to their commands,
inputs, negative controls and limitations rather than an inferred
fit or exclusion.
