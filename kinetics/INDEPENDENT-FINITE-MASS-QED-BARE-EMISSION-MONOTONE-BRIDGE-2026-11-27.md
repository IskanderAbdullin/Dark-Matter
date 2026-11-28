# Independent monotone emission integration resolves the soft real band

Extend only the **emission** sign of the
[signed moving-aperture bridge](INDEPENDENT-FINITE-MASS-QED-SIGNED-SUPPORT-BRIDGE-2026-11-27.md)
from the common signed interval to lab photon energies up to .00256.
This is the bare $P\leftrightarrow b+s+\gamma$ channel, with one
charge assignment, the same cold scalar momentum, and the same
frozen occupations and photon temperature 1. It is not an extension
of the monotonicity claim for the absorption sign.

At parent energy at least .18, $E_P/p$ is bounded by 1.042
for the archived parent mass. For positive photon energy the
minimum of the longitudinal photon derivative is $1-E_P/p$;
the scalar-angle derivative is bounded by $q E_P/p$. The
combined aperture equation is still strictly increasing on
[.18,6] for this emission band. The boundary equations are
inside the checked bracket for every sampled direction.
A bound on the available light-cone transfer excludes lower
energies; monotonicity excludes reentry above 6.
Shared-endpoint signed evaluations at photon energy .0002
coincide with the previous bridge. In the refactored solver
the signed defaults and the extended emission controls
reproduce the staged implementations bit-for-bit.

Separate bare gain and inverse loss are integrated directly
with their own photon Bose factors; there is no subtraction
of two huge cumulative cutoff totals. The six intervals
span 4e−5 to .00256. The largest two-rule relative difference
is 2.66e−13 (fine) and 7.18e−6 (third). For the first interval
[4e−5,8e−5], the unscaled gain/loss pairs are (65922.005,
19001.130) and (71187.538, 31061.389), respectively.
No hard anchor, old acceptance-window contribution, or
finite virtual term is included in these numbers.
In particular their sign and size do not determine the
matched thermal rate.

Replay from kinetics with pinned archived spectra:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_emission_monotone_bridge_reference.py
~~~
[Fine receipt](results/finite-mass-qed-bare-emission-monotone-fine-2026-11-27.json),
[third receipt](results/finite-mass-qed-bare-emission-monotone-third-2026-11-27.json).
