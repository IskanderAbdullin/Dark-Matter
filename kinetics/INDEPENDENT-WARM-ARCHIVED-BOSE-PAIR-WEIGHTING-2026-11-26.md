# Measure actual Bose weighting on five archived warm pair geometries

A maximum occupation need not be reached by the states that
a collision actually produces. Fold the stable fermion pair
annihilation amplitude twice on the same frozen geometry,
with and without scalar outgoing factors. Use q16 and q32
at their [matched 32-block early endpoint](INDEPENDENT-WARM-EARLY-SUPPORT-32-BLOCKS-2026-11-26.md),
and three later q16-style finite-grid archives from their
respective warm calculations. These are distinct checks,
not an interpolated or newly evolved trajectory.

For the vacuum fold, project both outgoing daughters onto
the same scalar energy shells and verify number and energy
in the projection. For each shell the exact within-shell
extrema of the reconstructed occupancy follow from its
entropy-variable endpoints. If `L≤f≤U` there, summing
vacuum daughter weights times the lower/upper extrema
encloses the linear Bose term. The quadratic term is
nonnegative and obeys
`fₐf_b ≤ f_peak(fₐ+f_b)/2`. This gives a separate check
of the directly folded increase on each positive geometry.
Unrecorded exterior scalar occupancies are set to zero,
not bounded by this procedure.

| Retained archive | T, MeV | Relative direct Bose increase in forward count | Peak-only upper increase |
| --- | ---: | ---: | ---: |
| q16, 32 blocks | 10.442 | 3.3220×10⁻⁵ | 7.5064×10⁻⁴ |
| q32, 32 blocks | 10.442 | 3.3162×10⁻⁵ | 7.5250×10⁻⁴ |
| warm b=3.5, 48 blocks | 3.434 | 1.7990×10⁻⁶ | 1.2211×10⁻³ |
| warm b=10, 42 blocks | 1.206 | 1.5047×10⁻⁹ | 1.4157×10⁻³ |
| warm b=30, 32 blocks | .414 | 0 in the represented fold | 1.4583×10⁻³ |

The two highest quadratures on nonzero boosts agree to
better than 0.04% relative in the enhancement. At b=30
the represented outgoing shell check also gives zero
enhancement; this is consistent with the earlier restricted
[accessible-band test](INDEPENDENT-WARM-PAIR-SUPPORT-CLOSURE-2026-11-26.md).
This says nothing about missing exterior scalar profiles,
intermediate transient peaks, the dynamics of other
processes or the complete signed coupled trajectory.
Do not turn these frozen phase-space ratios into an
occupation theorem for a redistributing dark sector.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_archived_bose_pair_weighting_reference.py \
  "$D"/warm-early-q{16,32}-time32-checkpoint-2026-11-26.json.gz \
  "$D"/warm-coupled-pair-b{35-48,10-42,30-32}-shape-2026-11-26.json
~~~
