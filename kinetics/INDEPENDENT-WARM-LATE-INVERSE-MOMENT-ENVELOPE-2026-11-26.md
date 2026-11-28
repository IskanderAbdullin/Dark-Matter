# Bound instantaneous mixed inverse gain without assigning a shape to known exterior receipts

Use the actual corrected q64-48 final checkpoint at 7.394049 MeV.
Keep the count and energy in every resolved shell, and the count,
birth mean and original momentum floor of every saved exterior
cohort. Transport those positive cohorts collisionlessly to the
anchor, then ask how large the **instantaneous** inverse-decay
source could be if all the recorded cohorts participated. The
result excludes any unrecorded population and makes no claim
about the distribution after reacting. Use the physical parent
cut .707792 GeV (energy .714821 GeV).

## Finite admitted support, moments only

For a scalar/light pair of energy intervals, first require
their sum to reach the parent cut. The exact two-body angular
condition also gives a lower energy of either leg as a convex
function of the other's energy. Its minimum on the remaining
partner interval is the unconstrained stationary point clipped
to that interval. This removes even the rest-endpoint singularity
from every allowed interval pair at the tested finite caps.
A numerical kinematic check kept all 492,847 sampled valid pairs.

For each leg, bound the eligible count above t using its
positive mean energy above the support floor, with the threshold
converted back to birth energy for an exterior cohort. In the
inverse-source angular measure the weights of the count are
bounded by 1/(p(t)t) and 1/p(t), respectively, for integrals
of occupation and of energy times occupation. Sum the pair
contributions, omitting heavy Pauli blocking conservatively.
For positive above-cut excess production use the merger
defect from the [convex tail audit](INDEPENDENT-WARM-CONVEX-HIGH-ENERGY-TAIL-2026-11-26.md),
bounded by the smaller of the two incoming energies and the cut.

| Assumed present momentum cap (GeV) | Bare inverse energy upper / retained ρH | Positive excess-energy upper / retained ρH |
| ---: | ---: | ---: |
| .70 | 1.3674e-5 | 5.3460e-7 |
| .85 | 3.3164e-5 | 6.7286e-7 |
| 1.05 | 3.4347e-5 | 6.7562e-7 |
| 1.50 | 3.5717e-5 | 6.7737e-7 |
| 2.00 | 3.6098e-5 | 6.7758e-7 |

These bounds do not assume the unknown shell shapes match the
reference entropy reconstruction. The moment-only two-core
inverse bound is 1.231e-29 in the same units, versus the exact
reference reconstruction 4.117e-32. Most of the all-sectors
upper comes from **one resolved and one exterior** leg; the
two-exterior allowance is only 1.650e-16 at the anchor.

## Unbounded exterior, with an explicit retained occupation ceiling

Drop the high momentum cap on the known exterior. In addition
to the moments, cap only each **resolved** occupation by the
extremum of its archived affine entropy profile. This controls
the near-rest resolved leg even if its exterior partner is
arbitrarily energetic. There is no pointwise occupation cap
assumed for the exterior. The corresponding upper values are
3.047e-5 and 6.689e-7. The rate normalization uses the
retained dark energy and H with the known exterior Jensen
gravity contribution; it does not add artificial active energy
to dilute the comparison.

As a cross-check and modest improvement, bound the mixed
source as a linear integral against each exterior birth-energy
distribution. Enclose its kernel on positive finite panels
and a final semi-infinite panel. For each trial nonnegative
slope choose the smallest intercept that majorizes every panel
(`alpha = max[K_upper − beta(E_birth,lower − floor)]`).
The known mean then bounds the rate by
`N*(alpha+beta*(mean−floor))`. Each trial is an upper, so
taking the minimum is safe without finding the optimal slope.
The finite-angle kernels were checked independently against
point evaluations up to 100,000 GeV. The refined bare bound is
2.974e-5; the excess bound is unchanged at 6.689e-7.

Neither calculation bounds the late evolution, recycling or
the lifetime of a frozen stress rate. In particular, the
retained ceiling is an anchor condition, not a proved
dynamical phase constraint. The bare source and positive
excess allowance also do not determine the **signed** flux.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/warm-corrected-anchor-signed-floor-2026-11-26.json
b=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
python -W error warm_late_mixed_inverse_moment_envelope_reference.py EOS.dat "$a" "$b"
python -W error warm_late_exterior_inverse_dual_reference.py EOS.dat "$a" "$b"
~~~

[Moment sector bounds](results/warm-late-inverse-moment-envelope-2026-11-26.json);
[semi-infinite exterior affine majorants](results/warm-late-exterior-inverse-dual-2026-11-26.json).
