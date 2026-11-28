# Why this contour envelope cannot be a stationary linewidth

The finite path calculation prepares both trajectories at the same
point, then lets the distance grow. Its early-time behavior is a useful
diagnostic of that preparation, not the short-time autocorrelation
of the production operator in a stationary bath.

On a fixed bounded soft momentum domain, the closed-current expansion
starts at order t²: the direct displacement cancels between all
segments, and the next term is the shrinking contour area. The included
spacelike and timelike symmetric contractions therefore have
χ(t)=O(t⁴). Matched t doubling at the smallest represented times gives
log₂ covariance ratios 3.9989 and 3.9956, using the normalized finite
covariance with its factored-out time restored. The covariance is
nonzero at longer time.

This precludes treating exp(−χ(t)) alone as the modulus of a nontrivial
stationary line with a nonnegative spectral weight. To see why, normalize
the putative nonnegative measure to unit area. The squared modulus of
its characteristic function is the characteristic function of the
difference of two independent frequencies from that measure. If its
departure from one is o(t²), Fatou's lemma applied to
(1−cos(frequency-difference times t))/t² says the second moment of
that difference is zero. Its only possible distribution is a point
mass at zero, so the original measure must be a single atom and its
modulus cannot have a later decay. Adding only a phase cannot cure
this obstruction. It is a statement about the isolated finite-cutoff
contour factor, not a restriction on the short-time behavior of the
true scalar spectral function.

An independent example with the time-white smooth spatial Gaussian
path noise explicitly gives a negative cosine-based density: minimum
−.04124 in angular-half-span units at sampled frequency 2.20. Extending
the integration time from 30 to 60 does not change this conclusion.
A well-defined stochastic expectation for a freshly initialized pair
does not automatically have a nonnegative stationary Fourier transform.

The next input has to be a matched neutral-operator correlation, with
its statistical and retarded components, the equilibrium detailed-balance
check where applicable, and appropriate open-system time evolution.
The finite-contour numbers are diagnostic constraints on such a
calculation, not a recipe for replacing the on-shell delta.

Reproduce with early_pole_wilson_stationary_line_gate_reference.py.
Raw data: results/early-pole-wilson-stationary-line-gate-2026-11-26.json.gz.
