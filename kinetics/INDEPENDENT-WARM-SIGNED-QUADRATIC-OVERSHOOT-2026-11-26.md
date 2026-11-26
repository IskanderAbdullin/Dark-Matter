# Test how far the signed bound can extend past one-half fermion occupation

The [time-resolved signed bound](INDEPENDENT-WARM-SIGNED-FERMION-TIME-COMPARISON-2026-11-26.md)
uses (f_1+f_2-1\le0) to drop the quadratic scalar Bose
term without imposing a scalar occupation ceiling. For a
uniform fermion per-spin cap ε above 1/2, its positive
part is at most (2ε-1). Multiply that excess by a
vacuum scalar-to-fermion pair rate with no final
blocking. This positively bounds the possible lost
fermion number from that piece of the signed operator;
it does not assume any monotonic Bose occupation.

Use detailed balance of the vacuum amplitudes and the
equal-mass energy-product weighting (E E'\ge s/4).
The same explicit finite invariant/angle cells and
a geometric high-invariant tail provide upper
coefficients for producing either fermion mass from
a scalar pair. The finer sums are
1.66204×10⁻⁸ and 2.09558×10⁻¹⁰ GeV⁻⁴, respectively,
when multiplying the squared physical scalar energy
density in the number-loss budget. Sum the two
coefficients; both use the same scalar reservoir.
Use the full dark energy upper for that scalar density
as a worst case. The resulting positive constant
term enters each time cell with (b^{-3}/H), alongside
the number-dependent losses.

At ε=1/2 the result reproduces the earlier comparison
exactly. For q32 the lower bound still exceeds .120
at ε=.50010 (about .396 in Ωh²), but this conservative
whole-energy treatment becomes inconclusive by
ε=.50014. A larger allowed cap is not evidence of
a physical high-occupation state or of depletion. This
is only a small extension of the sign range, and shows
where one needs either a much more selective phase-space
bound or a dynamical check on the cap. No scalar or
fermion spectrum is assumed from this comparison.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_signed_quadratic_overshoot_reference.py \
 "$D"/warm-isolated-fermion-survival-bound-2026-11-26.json \
 "$D"/warm-fermion-survival-with-leptons-2026-11-26.json \
 "$D"/warm-isolated-recycling-energy-envelope-2026-11-26.json \
 "$D"/warm-signed-fermion-time-comparison-2026-11-26.json
~~~
