# Halve the first observed q32 block without moving its endpoint

The [instrumented blocks](INDEPENDENT-WARM-OBSERVED-SUBOPERATOR-OCCUPATION-2026-11-26.md)
found a light occupancy maximum after the first pair half-kick
that was larger than both outer endpoints. Resume from the same
q32 source, replace the first coarse block by two of half the
logarithmic width, and keep the same final expansion and
10.18714 MeV temperature. Both runs keep their actual exterior
restart ledgers. Attach the same read-only observer.

| Observed q32 light peak | One coarse block | Two half-size blocks |
| --- | ---: | ---: |
| Shared entrance | 5.79567008×10⁻⁵ | 5.79567008×10⁻⁵ |
| Largest among all sampled interfaces | 5.80318938×10⁻⁵ | 5.79949004×10⁻⁵ |
| Shared-clock exit | 5.78616881×10⁻⁵ | 5.78568580×10⁻⁵ |
| Maximum excess over the two interval endpoints | 7.51930×10⁻⁸ | 3.81996×10⁻⁸ |

The within-interval excess falls to 50.8% of the
coarse value. In each of the two refined blocks, a
first pair half-kick still causes a small local rise
and the implicit update reduces it. This is a
clear discretization sensitivity, not a proof of
convergence or of a physical excursion in the
unsplit dynamics. In particular, endpoint samples
cannot certify a pointwise time envelope.

At the matched exit the normalized finite-shell
L1 differences of the scalar and stable number
moments are .0004731 and .0002900. The new upper
exterior creation is (4.06519×10^{-24}) versus
(4.06018×10^{-24}) GeV⁴, .123% smaller on the
refined interval. This is positive one-way receipt,
not reacting peripheral reuse. Old unknown tails,
the long warm evolution and the physical continuity
of the light occupation cap remain open.

For the fine tape, run the same observer twice
with --blocks 2 --stop 1.004987562112089, preserving
the exact source and the atomic state across passes.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_observed_suboperator_refinement_reference.py \
 "$D"/warm-observed-suboperator-q32-tape-2026-11-26.json.gz \
 "$D"/warm-observed-suboperator-q32-halved-first-block-2026-11-26.json.gz \
 "$D"/warm-light-occupation-next-q32-tape-2026-11-26.json.gz
~~~
