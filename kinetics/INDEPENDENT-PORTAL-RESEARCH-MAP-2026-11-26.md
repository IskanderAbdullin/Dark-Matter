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
   remain a separate check. The [matched 32-block follow-up](INDEPENDENT-WARM-EARLY-SUPPORT-32-BLOCKS-2026-11-26.md)
   closes all four time ladders on the same endpoint while
   documenting the separate fourth-moment sensitivity. The [q16 time ladder](INDEPENDENT-WARM-EARLY-TIME-REFINEMENT-2026-11-26.md)
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
The [separate isolated energy envelope](INDEPENDENT-WARM-ISOLATED-RECYCLING-ENVELOPE-2026-11-26.md)
does allow dark redistribution and repeated recycling in the
listed network after 10.442 MeV. With no new energy or reentry
and the nominal visible clock, it bounds the warm scalar density
to at most 2.35% if the total equals .120, after bounding the
virtual angular integral analytically. This still conditions the
bound on the starting energy and subsequent clock.
The [known-exterior extension](INDEPENDENT-WARM-KNOWN-EXTERIOR-RECYCLING-2026-11-26.md)
also allows all recorded old and new cohorts to participate
after the anchor, without assigning momenta to the old ones.
It leaves that fraction unchanged at the quoted precision.
A distinct [isolated fermion survival bound](INDEPENDENT-WARM-ISOLATED-FERMION-SURVIVAL-2026-11-26.md)
asks whether the *total* can fall to .120. Redistribution,
reverse gains and intermediate heavy fermions are allowed
but the outgoing scalar occupation must be capped by hand.
The [refined fermion floor with vacuum lepton loss](INDEPENDENT-WARM-FERMION-SURVIVAL-WITH-LEPTONS-2026-11-26.md)
uses exact angular ceilings and a finite energy budget for
high-invariant loss, retaining the explicit scalar cap.
A [finite-duration Bose test](INDEPENDENT-WARM-BOSE-WEIGHTED-DEPLETION-BUDGET-2026-11-26.md)
replaces that uniform cap by an expansion-weighted condition
and evaluates separate hypothetical burst lengths.
The [signed fermion-cap comparison](INDEPENDENT-WARM-SIGNED-PAIR-FERMION-CAP-2026-11-26.md)
does not assume a scalar cap. Instead a light and heavy
per-spin fermion cap controls the sign of the quadratic
Bose term. This alternative remains a conditional result.
Its [time-resolved lower comparison](INDEPENDENT-WARM-SIGNED-FERMION-TIME-COMPARISON-2026-11-26.md)
uses the number decrease within the analytic comparison,
giving a useful floor all the way to fermion cap one half,
still without claiming that cap dynamically.
A [small extension past one half](INDEPENDENT-WARM-SIGNED-QUADRATIC-OVERSHOOT-2026-11-26.md)
checks the previously nonpositive quadratic term. It shows
that the coarse positive allowance gives little extra range.
The [selected late charged-pion control](INDEPENDENT-WARM-SIGNED-FERMION-WITH-CHARGED-PIONS-2026-11-26.md)
adds that effective visible source, heat and reverse loss
to the signed floor. The momentum support matters for
the reverse bound; this is not a full hadronic sink.
The [heavy occupation barrier](INDEPENDENT-WARM-HEAVY-OCCUPATION-BARRIER-2026-11-26.md)
checks pointwise heavy creation under a prescribed light cap,
rather than inferring an occupation from its total rate.
It leaves the light cap and heavy elastic gains open.
A [direct light gain diagnostic](INDEPENDENT-WARM-LIGHT-DIRECT-PAIR-DOSE-2026-11-26.md)
locates the inconclusive term: conversion from the heavy
reservoir needs a much tighter energy fraction or source bound.
The separate [archived outgoing-pair Bose weights](INDEPENDENT-WARM-ARCHIVED-BOSE-PAIR-WEIGHTING-2026-11-26.md)
measure the boost on five actual retained shapes without
turning it into an evolving constraint.
The [late vacuum charged-lepton envelope](INDEPENDENT-WARM-LATE-LEPTON-RECYCLING-2026-11-26.md)
allows both visible pair injection and elastic heating in that
same period. At total .120 its upper scalar share is 2.44%.
It does not include late hadronic or new dark sources.
The [from-empty 50 MeV benchmark](INDEPENDENT-WARM-FROM-EMPTY-RECYCLING-ENVELOPE-2026-11-26.md)
is a separate numerical upper comparison from instantaneous
reheating. An energy-weighted vacuum kernel envelope and two
visible clocks give scalar ceilings below .060; no realistic
hadronic or additional scalar number source is included.
The [50 MeV selected pion-band and source headroom test](INDEPENDENT-WARM-PION-BAND-AND-SOURCE-HEADROOM-2026-11-26.md)
adds a prescribed vacuum KLOE exclusive window and a positive
unknown-source allowance. It does not bound unknown channels
or determine full thermal hadronic rates.
A [separate two-pion control](INDEPENDENT-WARM-PION-STATISTICS-AND-THRESHOLD-2026-11-26.md)
allows finite Bose enhancements, a capped missing threshold
strip and conditional pion elastic heating. The form-factor
caps there are model conditions, not data.
The [BABAR dispersive moment transfer](INDEPENDENT-WARM-BABAR-DISPERSIVE-GAP-2026-11-26.md)
adds a data-conditioned bound on just the 0.30–√.10 GeV
slice without inventing bins or duplicating KLOE coverage.
The [lower dispersive-theory threshold check](INDEPENDENT-WARM-DISPERSIVE-THEORY-PION-THRESHOLD-2026-11-26.md)
uses a separate conditional positive moment below .30 GeV,
not a substitute measured spectrum.
The [joint pion comparison](INDEPENDENT-WARM-JOINT-BABAR-PION-BUDGET-2026-11-26.md)
keeps the three charged bands disjoint and includes a separate
conditional neutral and elastic allowance. It is not full hadrons.
The [unmeasured-source budget by invariant band](INDEPENDENT-WARM-UNMEASURED-SPECTRAL-BUDGET-2026-11-26.md)
translates missing positive rates into sufficient conditional
caps and identifies where low-mass spectral coverage matters.
The [fixed-width coupling sensitivity](INDEPENDENT-WARM-VACUUM-COUPLING-SENSITIVITY-2026-11-26.md)
shows how changing the specified vector mixing and Yukawa
moves the 50 MeV sufficient upper boundary. It is not a fit.

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
