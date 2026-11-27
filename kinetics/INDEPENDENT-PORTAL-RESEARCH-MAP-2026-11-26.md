# Independent portal calculations: what each result does and does not test

An index to the independent numerical work in
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
   documenting the separate fourth-moment sensitivity.
   The [next short matched occupation stage](INDEPENDENT-WARM-LIGHT-OCCUPATION-NEXT-MATCHED-STAGE-2026-11-26.md)
   continues both by a further 2% in expansion, auditing
   the exact saved occupation extrema at every new block.
   A [further instrumented warm continuation](INDEPENDENT-WARM-OBSERVED-SUBOPERATOR-OCCUPATION-2026-11-26.md)
   shows small increases within blocks despite descending
   saved endpoint light maxima.
   The [matched q32 half-block check](INDEPENDENT-WARM-OBSERVED-SUBOPERATOR-REFINEMENT-2026-11-26.md)
   reduces the observed internal increase by nearly half
   without claiming a continuous occupation bound. An independent
   [two-moment feasibility audit](INDEPENDENT-WARM-SHELL-OCCUPATION-FEASIBILITY-2026-11-26.md)
   proves that the saved light shell moments fit both very dilute
   profiles and profiles reaching Pauli saturation; it is not
   legitimate to deduce a dynamical cap from those moments. A [conditional shell regularity test](INDEPENDENT-WARM-SHELL-LIPSCHITZ-CAP-2026-11-26.md)
   quantifies the additional physical-momentum slope bound
   that would suffice at the saved times. A [frozen pair-fold counterexample](INDEPENDENT-WARM-SHELL-SATURATED-PAIR-FOLD-2026-11-26.md)
   compares different profiles with the same shell moments; a change
   in the largest-number shell moves the net pair rate by about 1%. The [near-rest inverse-decay check](INDEPENDENT-WARM-NEAR-REST-LIGHT-INVERSE-DECAY-2026-11-26.md)
   localizes the innermost light response to a far-tail scalar window;
   cropping only that support flips the frozen q32 rate sign. A [selected visible injection bound](INDEPENDENT-WARM-COLD-LIGHT-VISIBLE-MARGIN-2026-11-26.md)
   leaves the sampled full-support decay-only sign negative,
   without certifying the sign of all coupled operators. The [tagged dark-pair follow-up](INDEPENDENT-WARM-COLD-LIGHT-PARENT-DARK-PAIR-2026-11-26.md)
   resolves tiny represented positive gains into the relevant heavy
   modes, on those fixed shapes only. A [shell-number envelope](INDEPENDENT-WARM-COLD-PARENT-RETAINED-PAIR-ENVELOPE-2026-11-26.md)
   bounds the extra retained dark pair gains without picking
   their within-shell shapes; the baseline rate stays frozen. The [tagged light scalar-pair fold](INDEPENDENT-WARM-COLD-LIGHT-SCALAR-PAIR-2026-11-26.md)
   resolves a second negative cold contribution, much less
   sensitive to the same scalar support crop. A [below-pole light elastic fold](INDEPENDENT-WARM-COLD-LIGHT-BELOW-POLE-ELASTIC-2026-11-26.md)
   finds a separate local sink of comparable size while omitting
   the pole and the above-cut region. A [whole-shell elastic projection audit](INDEPENDENT-WARM-COLD-ELASTIC-PROJECTION-AUDIT-2026-11-26.md)
   resolves the remaining first-cell signed-rate sensitivity
   and compares the native and independent quadratures. The [native prompt projection comparison](INDEPENDENT-WARM-COLD-PROMPT-SUPPORT-GRID-2026-11-26.md)
   quantifies the disconnected q16 cold shells, and independently
   reproduces the q32 first-cell decay rate. The [matched open prompt boundary check](INDEPENDENT-WARM-POST-CHECKPOINT-OPEN-PROMPT-FLUX-2026-11-26.md)
   finds opposite scales of frozen spill on the two supports;
   cropping the same q32 core reproduces the q16 large rate. The [explicit reacting-exterior transient](INDEPENDENT-WARM-PROMPT-REACTING-EXTERIOR-TRANSIENT-2026-11-26.md)
   evolves both species for a very short fixed-clock test; its actual
   outside counts are much smaller than a one-way linear extrapolation. A [matched tail-shape audit](INDEPENDENT-WARM-PROMPT-REENTRY-SHAPE-AUDIT-2026-11-26.md)
   finds that the newly populated tail is close in total moments
   and normalized shell shape to the actual q32 checkpoint. The [discrete gross-crossing ledger](INDEPENDENT-WARM-PROMPT-GROSS-CROSSING-2026-11-26.md)
   explicitly measures the inward and outward components; their
   difference recovers the tiny resident number. A [short expanding prompt-only control](INDEPENDENT-WARM-PROMPT-REENTRY-EXPANDING-2026-11-26.md)
   extends the matched support to q64 and tracks expansion work
   as well as the represented tail. The [coupled support comparison](INDEPENDENT-WARM-MATCHED-EXPLICIT-SUPPORT-COUPLED-2026-11-26.md)
   repeats the test on a slightly longer matched expansion with
   selected elastic and signed pair reactions and raw checkpoints. A [known-cohort activation control](INDEPENDENT-WARM-REACTIVATE-KNOWN-COHORTS-2026-11-26.md)
   allows stored moment packets to react without counting
   them twice in gravity, for one stated re-entry prescription. A [four-level q32 time ladder](INDEPENDENT-WARM-Q32-MATCHED-TIME-LADDER-2026-11-26.md)
   observes near-first-order temporal convergence through
   thirty-two actual blocks at the same physical endpoint. A [moment extrapolation diagnostic](INDEPENDENT-WARM-Q32-TIME-EXTRAPOLATION-2026-11-26.md)
   cancels much of the leading difference while preserving
   reconstruction admissibility on this ladder only. A [borrowed-guess replay](INDEPENDENT-WARM-Q64-TRANSPLANTED-NEWTON-2026-11-26.md)
   reproduces a wider implicit block by using the archived
   narrower shape only to initialize the nonlinear solve. A [three-interval q64 continuation](INDEPENDENT-WARM-EXTENDED-COUPLED-SUPPORT-2026-11-26.md)
   matches the already archived q32 time grid through 10.238 MeV. A [guided q64 continuation](INDEPENDENT-WARM-Q64-GUIDED-SERIAL-2026-11-26.md)
   advances two more accepted wide blocks through 10.137 MeV,
   using q32 only as initial guesses. A [q64 time refinement](INDEPENDENT-WARM-Q64-GUIDED-TIME-REFINEMENT-2026-11-26.md)
   keeps wider physical equations fixed while comparing two,
   four and eight blocks to the same temperature. A [history-only q64 guess](INDEPENDENT-WARM-Q64-HISTORY-PREDICTOR-2026-11-26.md)
   reproduces a wider continuation block without using the next
   narrower shape in that solve. A [history-only q64 series](INDEPENDENT-WARM-Q64-HISTORY-SERIAL-2026-11-26.md)
   advances seven further blocks and verifies the atomic
   restart and number ledger through 9.748 MeV. An independent q32
   continuation to the same endpoint checks the shared shells. A [matched later continuation](INDEPENDENT-WARM-HISTORY-LOWER-MATCHED-2026-11-26.md)
   advances both supports independently through another sixteen
   accepted blocks to 9.027 MeV. The matched 32-block test
   separates time sensitivity from the much smaller support change. A [matched cold recoil crossing](INDEPENDENT-WARM-RECOIL-CROSSING-MATCHED-2026-11-26.md)
   then checks which local prompt channel vanishes on q32 and
   which remains on q64 around 8.24 MeV. A [finite-momentum tagged pair audit](INDEPENDENT-WARM-COLD-LIGHT-FINITE-PAIR-2026-11-26.md)
   finds that the native reverse quadrature materially undercounts gain into the first cold cell. The [conservative cold-event partition](INDEPENDENT-WARM-COLD-PAIR-PARTITION-2026-11-26.md)
   removes a one-step native cold spike without shifting the full spectrum materially in that block. A [six- and nine-shell lift](INDEPENDENT-WARM-COLD-PAIR-MULTISHELL-2026-11-26.md)
   checks adjacent cold cells without changing the forward operator or the remaining inverse phase space. A [24-step corrected q32 continuation](INDEPENDENT-WARM-MULTISHELL-CORRECTED-TRAJECTORY-2026-11-26.md)
   tests the same change over the whole matched 8.133 to 7.394 MeV interval.
   The [native recoil-gate time/support cross](INDEPENDENT-WARM-NATIVE-RECOIL-CROSS-AND-GATE-LAW-2026-11-26.md)
   checks the coarser trajectory and derives the local quadratic shutoff
   from the narrowing light-momentum/recoil-energy region.
   A [corrected 12/24-step comparison](INDEPENDENT-WARM-CORRECTED-Q32-TIME-LADDER-2026-11-26.md)
   measures the much larger whole-state time sensitivity and checks that the cold-cell decline survives it. A [thermal balance and endpoint budget check](INDEPENDENT-WARM-LIFTED-PAIR-BALANCE-2026-11-26.md)
   independently tests gain/loss cancellation in narrow cells and resolves the remaining elastic contribution.
   A [common-coordinate pair kernel](INDEPENDENT-WARM-PAIR-COMMON-CHART-2026-11-26.md)
   tests both directions on one phase-space measure, replays the q32 spike at two orders,
   and solves one wider q64 replay. A [number/heat Onsager control](INDEPENDENT-WARM-COMMON-PAIR-NUMBER-HEAT-2026-11-26.md)
   tests full quantum pair transfer between two near-equilibrium baths. A [frozen scalar-pair sign bound](INDEPENDENT-WARM-FROZEN-PAIR-BARRIER-2026-11-26.md)
   separates represented blocking from the suppressed positive vacuum-partner channel, without establishing a live fermion cap.
   The [completed native q32/q64 baseline](INDEPENDENT-WARM-BEYOND-RECOIL-NATIVE-2026-11-26.md)
   separates the recoil cutoff from the native common pair-projection spike.
   [Entropy perturbation envelopes](INDEPENDENT-WARM-PAIR-ENTROPY-ROBUSTNESS-2026-11-26.md) test which
   perturbations preserve the selected pair sign without treating observed numerical discrepancies as true errors.
   [Lower frozen sign thresholds and thermal exchange modes](INDEPENDENT-WARM-PAIR-THRESHOLD-AND-MODES-2026-11-26.md)
   distinguish the local sign argument from two-reservoir number/heat relaxation.
   [Same-moment thermal and exterior scalar tests](INDEPENDENT-WARM-PAIR-PHASE-UNCERTAINTY-2026-11-26.md)
   quantify sampled shape sensitivity and give conditional nonvacuum exterior bounds.
   The [complete corrected time/support ladder](INDEPENDENT-WARM-CORRECTED-TIME-AND-SUPPORT-2026-11-26.md)
   compares 12, 24 and 48 q32 steps plus 24 wider q64 steps. It
   isolates the time error from the tested outer-support change;
   the cold decline survives both controls.
   The [completed four-run time/support cross](INDEPENDENT-WARM-CORRECTED-FOUR-RUN-TIME-SUPPORT-2026-11-26.md)
   adds q64 at 48 steps and checks the signed mixed effect separately.
   The [late complete common-chart replay](INDEPENDENT-WARM-Q64-LATE-WHOLE-PAIR-REPLAY-2026-11-26.md)
   preserves the cold decline at another interval, but its fallback
   prevents a separate attribution of all differences to quadrature.
   A [locked two-substep follow-up](INDEPENDENT-WARM-Q64-LATE-COMMON-LOCKED-CLOCK-2026-11-26.md)
   separates pair-order variation from the stiff inner clock.
   A [driven heavy-lag test inside two fine intervals](INDEPENDENT-WARM-Q64-HEAVY-LAG-ENTROPY-SURROGATE-2026-11-26.md)
   recalculates the parent source on an explicit daughter entropy interpolation.
   A [realizable first-order cancellation diagnostic](INDEPENDENT-WARM-CORRECTED-RECONSTRUCTED-TIME-DIAGNOSTIC-2026-11-26.md)
   reconstructs admissible Fermi/Bose phases and recomputes pressure,
   rather than advancing or directly interpolating those phases.
   The [completed frozen q64 pair test](INDEPENDENT-WARM-Q64-COMPLETE-FROZEN-PAIR-CONTROL-2026-11-26.md)
   checks all 25 stored states and the 24 prescribed interpolation brackets.
   The [actual four-interval split budget](INDEPENDENT-WARM-Q64-ACTUAL-COLD-OPERATOR-BUDGET-2026-11-26.md)
   compares the accepted pair, prompt and elastic contributions with the
   unchanged comoving cold average across collisionless redshift.
   An [ordered collision-only common thermal tangent](INDEPENDENT-WARM-Q64-COMMON-TANGENT-TRANSPORT-2026-11-26.md)
   follows the two prescribed transfer modes across five q64 anchors
   without turning the fitted reservoirs into an actual kinetic history.
   A [common equilibrium at the archive densities](INDEPENDENT-WARM-PAIR-COMMON-MOMENT-EQUILIBRIUM-2026-11-26.md)
   preserves aggregate number and energy and shows sub-Hubble slow pair relaxation in the cold reservoir comparator.
   The [q16 time ladder](INDEPENDENT-WARM-EARLY-TIME-REFINEMENT-2026-11-26.md)
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
A [piecewise radiation-clock comparison](INDEPENDENT-WARM-SIGNED-FERMION-PIECEWISE-CLOCK-2026-11-26.md)
updates the conditional survival reference at 9.748 and 9.027 MeV; it does
not establish the fermion cap. The [matched corrected-anchor sensitivity](INDEPENDENT-WARM-CORRECTED-SIGNED-FLOOR-SENSITIVITY-2026-11-26.md)
shows that the tested time/support differences barely move the same
conditional lower floor at 7.394 MeV. The [heavy occupation barrier](INDEPENDENT-WARM-HEAVY-OCCUPATION-BARRIER-2026-11-26.md)
checks pointwise heavy creation under a prescribed light cap,
rather than inferring an occupation from its total rate.
It leaves the light cap and heavy elastic gains open.
A [direct light gain diagnostic](INDEPENDENT-WARM-LIGHT-DIRECT-PAIR-DOSE-2026-11-26.md)
locates the inconclusive term: conversion from the heavy
reservoir needs a much tighter energy fraction or source bound.
The [aggregate finite-band heavy energy envelope](INDEPENDENT-WARM-HEAVY-COMOVING-BAND-ENERGY-2026-11-26.md)
bounds that reservoir below an explicit comoving cut,
without treating the higher band as absent.
The [late partial direct-gain recalibration](INDEPENDENT-WARM-LATE-DIRECT-LIGHT-GAIN-HEADROOM-2026-11-26.md)
updates the headroom and transports the old finite bands without
resetting them to empty or asserting the global light cap.
A [frozen high-inverse shape stress](INDEPENDENT-WARM-HIGH-INVERSE-FROZEN-PRESSURE-TEST-2026-11-26.md)
shows why total energy and a light cap do not force a uniformly tiny
source-to-decay reservoir without information about momentum location.
The [convex excess-energy ledger](INDEPENDENT-WARM-CONVEX-HIGH-ENERGY-TAIL-2026-11-26.md)
separates high-parent circulation from actual new above-cut energy and
bounds the latter for the known ballistic receipts without assuming a shape.
The [late fixed-moment reactivation stress](INDEPENDENT-WARM-LATE-RECORDED-TAIL-REACTIVATION-2026-11-26.md)
shows that even known stored moment bounds do not determine the initial
reacting rate or its sign when the unknown phase shape is changed.
The [late mixed-source envelope](INDEPENDENT-WARM-LATE-INVERSE-MOMENT-ENVELOPE-2026-11-26.md)
provides an instantaneous upper on known-cohort contributions with
finite support, or with unbounded exterior under a retained phase ceiling.
The [cut-aligned prompt control](INDEPENDENT-WARM-CUT-ALIGNED-PROMPT-VACUUM-CONTROL-2026-11-26.md)
keeps the physical vacuum in an accepted root and shows why an
initial retained phase ceiling need not persist after reactivation.
The [aligned .68 time refinement](INDEPENDENT-WARM-LATE-TAIL-TIME-REFINEMENT-2026-11-26.md)
tracks the transient decay and distinguishes full energy above the cut
from the smaller convex excess.
A [dynamical threshold-collar split](INDEPENDENT-WARM-LATE-TURNING-LOCAL-COLLAR-2026-11-26.md)
checks local momentum sensitivity at the same six turning times.
A [frozen pad-opening audit](INDEPENDENT-WARM-LATE-PROMPT-PAD-OPENING-2026-11-26.md)
checks direct outgoing support without claiming dynamical re-entry.
The [aligned .80 control](INDEPENDENT-WARM-LATE-ALIGNED-HARD080-RELAXATION-2026-11-26.md)
starts above the cut and has negative transfer at every saved time.
A [species-resolved comparison](INDEPENDENT-WARM-LATE-TAIL-SPECIES-EXCHANGE-2026-11-26.md)
tracks daughter currents. A [near-cut width evolution](INDEPENDENT-WARM-LATE-NEAR-CUT-SHAPES-2026-11-26.md)
finds opposite excess changes for matched initial moments. The [local expanding comparison](INDEPENDENT-WARM-LATE-COMPETING-FLUX-2026-11-26.md)
checks selected pair and elastic contributions for .68, with an
[already above-cut .80 control](INDEPENDENT-WARM-LATE-HARD080-COMPETING-FLUX-2026-11-26.md). An [open elastic receipt check](INDEPENDENT-WARM-LATE-OPEN-ELASTIC-RECEIPT-2026-11-26.md)
measures the separate outgoing-only channel. A [matched vacuum pad](INDEPENDENT-WARM-LATE-MATCHED-VACUUM-PAD-2026-11-26.md)
checks the bulk effect.
An [implicit finite-parent replay](INDEPENDENT-WARM-LATE-FINITE-PARENT-TURN-2026-11-26.md)
checks startup and confirms the late turn on the same discrete clocks. A [comoving parent chart](INDEPENDENT-WARM-COMOVING-PARENT-CHART-2026-11-26.md)
audit verifies redshift work and collision geometry. A [matched expanding test](INDEPENDENT-WARM-FINITE-PARENT-COMOVING-EXPANSION-2026-11-26.md)
combines finite lifetime, expansion and selected competing dark reactions. The [coupled clock audit](INDEPENDENT-WARM-COUPLED-CLOCK-AND-SPECIES-BUDGET-2026-11-26.md)
separates time resolution from the daughter species energy exchange. A [pair-versus-elastic comparison](INDEPENDENT-WARM-COUPLED-PAIR-ELASTIC-ATTRIBUTION-2026-11-26.md)
separates their opposite scalar energy increments. A [dynamic outgoing pad](INDEPENDENT-WARM-COMOVING-PROMPT-OUTGOING-2026-11-26.md)
resolves new daughter occupation without an incoming floor. A [two-stage stiff clock](INDEPENDENT-WARM-FINITE-PARENT-SDIRK-CLOCK-2026-11-26.md)
reduces .68 sensitivity but rejects a tested .80 stage. The [durable restart check](INDEPENDENT-WARM-DURABLE-COMOVING-RESTART-2026-11-26.md)
verifies byte identity for both accepted clocks.
The [conditional high-band closure](INDEPENDENT-WARM-HEAVY-HIGH-BAND-FINITE-DAUGHTER-SUPPORT-2026-11-26.md)
covers that part only with a separate persistent daughter
support hypothesis; visible thermal tails are still bounded. A [finite ancestry extension](INDEPENDENT-WARM-HEAVY-PERIPHERAL-ANCESTRY-2026-11-26.md)
excludes that band from one internally made exterior input
interacting with one retained input; two exterior inputs
or repeated re-entry are not bounded by that q32 test.
The [observed one-way exterior creation ledger](INDEPENDENT-WARM-POSITIVE-EXTERIOR-RECEIPT-2026-11-26.md)
measures new outer cohort energy on the saved 2% stage.
It is not gross reacting use if re-entry is switched on. The [eligible receipt tail](INDEPENDENT-WARM-HEAVY-ELIGIBLE-EXTERIOR-RECEIPT-2026-11-26.md)
uses the stored support and birth-energy mean to bound the
one-use high-momentum subset, even for two-exterior events.
That subset still is not a measured gross reacting budget. A [high-parent prompt probe](INDEPENDENT-WARM-PROMPT-HIGH-PARENT-LEDGER-2026-11-26.md)
measures the bare inverse-decay contribution on a separate short
explicitly reacting q64 control, not the full high-heavy budget. A [known-cohort shape-independent initial upper](INDEPENDENT-WARM-HIGH-INVERSE-KNOWN-COHORT-SHAPE-BOUND-2026-11-26.md)
retains the finite admitted support but does not pick an unknown
within-cohort distribution. A [resolved-shell moment envelope](INDEPENDENT-WARM-HIGH-INVERSE-RESOLVED-SHELL-ENVELOPE-2026-11-26.md)
independently bounds the later retained bare high inverse rates
without choosing their within-shell profiles. A [direct-occupation hard-tail stress](INDEPENDENT-WARM-DIRECT-PROMPT-STRESS-2026-11-26.md)
evolves two moment-matched hard packets briefly and separates
robust strong relaxation from unresolved terminal precision. The [isolated exposure time ladder](INDEPENDENT-WARM-DIRECT-EXPOSURE-TIME-LADDER-2026-11-26.md)
separates path error from integration samples and does not claim a comoving dose. A [completed selected coupled stress](INDEPENDENT-WARM-COUPLED-COHORT-STRESS-2026-11-26.md)
checks a hard .58 GeV realization on the same two-step grid as a vacuum-padded control
and resolves relaxation to the selected-network terminal rate in this case.
The [matched .52/.58 GeV comparison](INDEPENDENT-WARM-COUPLED-COHORT-HARDNESS-2026-11-26.md)
separates terminal agreement of the retained spectra from the orders-of-magnitude
sensitivity of the tiny new exterior receipt counts.
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
