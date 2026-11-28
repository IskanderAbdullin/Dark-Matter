# A cold threshold can move without any reactions

Let the two corrected 7.394 MeV archives coast in the same
[fixed reference background](INDEPENDENT-ARCHIVED-COLLISIONLESS-STREAMING-2026-11-26.md).
Hold all comoving number and momentum fixed, but account for
photon heating from the visible entropy table. The initial
momentum corresponding to p/T<.1 at reference equality is
1.397655 times the anchor threshold. In either archive about
1.8755e−4 of represented scalars meet the cut at the anchor,
and 4.9758e−4 meet it at equality, 2.65 times as many even
though no particles are created. The mean distances traveled
to equality by those two selected scalar subsets are .533
and .712 Mpc respectively. The full represented scalar mean
is 10.418 Mpc, while the stable species mean is .213 kpc.

The represented scalar contributes only 2.6324e−6 of the
two-species rest mass in these frozen archives. Known exterior
numbers alter this fraction by less than 2.3e−11 relative
without choosing an exterior momentum spectrum. A separate
free-propagation monopole uses j0(k chi) for an initially
common fractional perturbation and fixed rest-mass weighting.
Compared with stable-only, adding scalars changes that
dimensionless monopole by at most 2.64e−6 on the eight sampled
modes .003–100 Mpc⁻¹, and by at most twice the small scalar
fraction at any k in this particular control. This uses
no gravitationally driven perturbation, baryon response,
or growing mode.

Refine angular phase in the radial quadrature without
changing the input occupations, and refine the expansion
quadrature separately. A first unsplit result showed why
number/mean distance convergence cannot establish convergence
of the oscillatory large-k scalar monopole; the phase-split
results agree on both full and cold-band monopoles. These
diagnostics distinguish temperature labeling, represented
abundance and kinematics, and are not a self-consistent
matter transfer function or an observational exclusion.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_late_cold_scalar_ballistic_impact_reference.py EOS.dat \
 results/warm-corrected-q32-time48-step48-2026-11-26.json.gz \
 results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
~~~

[Ballistic observables and phase-split convergence](results/warm-late-cold-scalar-ballistic-impact-2026-11-26.json).
