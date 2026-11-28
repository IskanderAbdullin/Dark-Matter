# Positive forward depletion, and a local early inverse-reaction check

This augments the [first-order recoil tally](INDEPENDENT-FINITE-COOLING-PAIR-RECOIL-2026-11-26.md).
It separates a mathematically bounded **forward-only** model from a
probe of whether dropping the inverse channel is justified.
The cooling clock, matching condition and external source are prescribed.

## Bounds for the forward-only split operator

Project every stable source bin onto its two fixed comoving momentum
endpoints, preserving number and its squared-momentum integral.
Within each reaction substep the nonnegative kernel \(K\) includes the
full angular threshold and the visible expansion clock. If \(Y^0\)
is the passive upper number vector, \(Y\) the depleting vector, and
\(L\) the lower vector at substep entry, then
\[
L_i e^{-s (K Y^0)_i}\le Y_i(s)\le Y^0_i,\qquad 0\le s\le1.
\]
Propagate the linear absorbing lower vector across source injections.
For positive quadratic emission, integrate the lower product exactly
using \((1-e^{-B})/B\), \(B=(K Y^0)_i+(K Y^0)_j\). The same bound
works with the nonnegative outgoing squared-momentum kernel. These
bounds apply to the exact **split finite-node** dynamics, not to the
unprojected continuum or to the midpoint approximation itself.

The positive midpoint trajectory stays between the bounds for both
scalar number and its momentum moment. From the refined projection,
stable \(\Omega h^2\) decreases from 0.120 to approximately 0.119930
(instantaneous) and 0.119882 (finite). Relative to the passive
rate with the **same** projection, feedback lowers new scalar counts
by about 0.063% and 0.16%, respectively. Differences from the
alternative within-bin source representation can exceed those small
feedback corrections. Conservation of the outgoing scalar number
against stable loss and positivity are checked independently. This
does not bound stimulated or inverse conversion.

## A necessary early signed-rate check

At labeled epochs, reconstruct the positive prompt spectra independently
with two log widths and evaluate lab-angle forward and inverse rates
including the occupation factors. At 25 MeV in the finite cooling
history the inverse/forward ratio is about 0.89–0.90, even though
the scalar occupation remains small. By 15 MeV the ratio is about
0.019, by 9 MeV roughly \(4\times10^{-6}\), and by radiation matching
roughly \(2\)–\(3\times10^{-10}\). The reverse channel is kinematically
open at the earlier hot epochs on these shapes; its absence in the
late-anchor probe cannot establish its early absence. The small
effect of forward stable depletion does **not** justify dropping
this inverse process during the whole production history.

These are local ratios on passive prompt-only shapes. The newly
produced recoil distribution is not fed back into either rate.
No signed evolving relic count or physical transfer follows from
this test alone.

As a separate diagnostic, integrate the spontaneous **signed** vacuum
kernels on the same evolving but forced prompt spectra. The sampled
local inverse/forward ratio can exceed one near the hot start.
Nevertheless the refined pre-radiation integrated forward count
contributes Ωh²=3.535×10⁻¹⁰ and the reverse removes 4.320×10⁻¹²,
about 1.22% of that pre-radiation term. The post-matching forward
term is another 1.513×10⁻¹⁰ with negligible prompt inverse in this
controlled fold. Total absorption is about 0.86% of the passive
forward total. Timing and redshift weights are essential: a large
ratio in one short early interval is not a large integrated effect.
An independent lab-angle normalization agrees to 3.85×10⁻¹¹
relative on the shared test grid. The fold still does not include
recoil feedback, depletion or evolving inverse occupations.

## Signed momentum correction on the same passive source

In each inverse event subtract the squared physical momenta of
the two incoming scalars. Although the inverse threshold requires
a large instantaneous physical momentum, the early events redshift
substantially before the 1 MeV anchor. Using the same signed
vacuum quadrature for both directions, the refined early absorbed
population has RMS \(p/T=0.432\) at the anchor; pre-radiation
emission has RMS 2.872 and later emission 7.140.
Relative to the sum of prompt and spontaneous forward populations,
the reverse subtracts 0.158% of counts but only 0.00705% of the
squared-momentum yield. The normalized RMS moves slightly upward,
2.0448 → 2.0464. This is a net integrated moment on a forced
source, not a reconstructed positive signed PSD or a coupled result.
## Projected passive removal in the tracer

The inverse number and squared momentum can also be deposited into
the same positive q bins as the prompt distribution before subtraction.
Their sums reproduce the integrated inverse rates and moments.
Using the projected removal as a separate linear LOS component
on the same fixed metric, the inferred correction to the forward-only
tracer at z=20 is −0.0001241 to −0.0001280 at k=4 and −0.0001982
to −0.0002044 at k=10 h/Mpc in the two original controls.
This effect is much smaller than the fluid-closure artifact.

The removed subpopulation's CDF can differ by up to 0.091 between
these grids. A third control using 104 incoming bins and 84
time panels changes the inverse CDF by 0.0345 relative to the
second control and gives combined shifts of −0.0001264 and
−0.0002016 at k=4 and 10. The change to the combined tracer is stable at about
6×10⁻⁶ absolute at k=10, not proof of a converged absorption
shape. The reconstructed joined-minus-removal spectrum was tested
at its knot union; this is not a continuous positivity proof.
There is still no evolving pair feedback in this subtraction.

## Reproduction

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_forward_depletion_reference.py \
    EOS.dat results/finite-cooling-reheat-portal-2026-11-26.json
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_forward_depletion_reference.py \
    EOS.dat results/finite-cooling-reheat-portal-2026-11-26.json --high
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_source_signed_probe_reference.py \
    EOS.dat results/finite-cooling-reheat-portal-2026-11-26.json
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_signed_passive_fold_reference.py \
    EOS.dat results/finite-cooling-reheat-portal-2026-11-26.json
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_signed_passive_fold_reference.py \
    EOS.dat results/finite-cooling-reheat-portal-2026-11-26.json --recoil-moments
~~~

For the projected removal add --absorption-grid followed by the
published finite-cooling-reheat-distribution result to the signed
fold with --recoil-moments. Add --refine-absorption for the third control. Then use the saved compact fixed metric:

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_inverse_recoil_los_reference.py \
  EOS.dat results/finite-cooling-passive-absorption-shape-2026-11-26.json \
  MIXED.dat results/finite-cooling-los-off-portable-2026-11-26.json \
  results/los-metric-off/*.npz
~~~
