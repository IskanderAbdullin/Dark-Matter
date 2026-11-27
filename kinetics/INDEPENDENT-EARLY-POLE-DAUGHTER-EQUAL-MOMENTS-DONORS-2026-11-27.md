# Equal daughter bare moments with different compensation supports

Keep the two aperture-local daughter cells exactly as substituted;
restore full retained daughter number and bare energy to their
native values by retuning different outside subsets. The same
two-parameter affine FD correction is used with the same 0.2
coefficient gate. Other cells, hard geometry and cold seed are
unchanged. Compare all outside cells, the lower side only, the
upper side only, and both sides after excluding the immediately
adjacent cell on each side.

| Allowed compensation support | Initial local-minus-native cold-rate shift | Orthogonal coefficient | Outcome |
|---|---:|---:|---|
| All outside | -0.0155110728714 | +0.05337 | Accepted |
| Below the local aperture | -0.0155110728199 | +0.05339 | Accepted |
| Above the local aperture | Not reported | Not credited | Failed moment and amplitude gates |
| Both sides, immediate neighbors excluded | -0.0155110599944 | +0.06825 | Accepted |

Accepted number and bare-energy residuals are at most 3.47e-18.
Excluding the two immediate neighbors increases the largest
compensator moment-component shift from about 5.18e-5 to 6.94e-5
in the same coordinates, but the initial cold rate remains close.
This static equality does not determine the finite evolution,
nor does the rejected upper-side fit prove no upper-supported
shape can match the target. It rejects this affine recipe under
the stated gates. A finite comparison of the two-sided supports
tests preparation shape while holding both bare moments fixed.

Run `early_pole_daughter_equal_moments_donor_support.py` against
the matched fine projected raw archive. The default preparation
path without a donor subset replays the previous static receipt
unchanged.
