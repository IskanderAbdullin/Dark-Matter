# A scheme-defined vacuum normalization with unresolved photons

The [relative vacuum vertex](INDEPENDENT-FINITE-MASS-QED-VACUUM-SCALAR-VERTEX-SHAPE-2026-11-27.md)
cancels the scalar-current normalization and the two external
residues. Restore them here with on-shell fermion masses/legs
and an MS-bar Yukawa coupling at mu=M. The calculation is
for a single Abelian gauge interaction in vacuum. Include
unresolved real radiation only below the stated energy cap
in the parent rest frame, and one virtual photon correction.

Obtain each external residue from the vector and scalar parts of
the one-loop fermion self-energy. With x the internal fermion
Feynman parameter and lambda the auxiliary photon mass, its
on-shell denominator is x² m²+(1−x) lambda². The d=4−2ε
coefficient in the vector numerator gives a finite +1
after UV subtraction. The virtual scalar L² coefficient
similarly supplies a finite constant beyond the relative
vertex. The UV poles from vertex plus external legs sum to
3 in units g²/(16 pi²), removed by the scalar MS-bar
counterterm. Before the real part, the remaining common
soft-log coefficient in the normalized rate is −B, the
opposite of the leading real angular radiator.

The residue integral agrees with an independent 66-digit
self-energy derivative. At equal fermion masses and zero
transfer the amplitude coefficient is −2 at mu=m, as
required by differentiating the one-loop Abelian pole mass
with respect to the running scalar mass insertion (the
one-loop Abelian color factor limit of the
[pole-to-running mass relation](https://arxiv.org/abs/hep-ph/9911434)).
The one-loop finite amplitude changes with scale by
6 ln(mu/M), giving the explicit scale dependence
needed to cancel the one-loop running of the Yukawa coupling.
The two unequal-mass finite window results approach the
regulator-free answer on all three photon-mass ladders.

At the two archived scalar poles and rest photon cap
omega/M=.01 the absolute vacuum coefficient is about −11.17,
or −7.07e−6 relative to Born at g=.005. The coefficient varies
from about −19.53 to −5.63 as the cap changes from .001 to .05.
This cut dependence is expected. The regulator used in the
checks applies to the vertex, external legs, and eikonal
real part; the finite real recoil remainder is evaluated
in its massless limit. None of these results includes
resolved hard photons, a thermal vertex, medium masses
or a calibrated cosmological production rate.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_absolute_vacuum_window_reference.py
```
[Raw receipt](results/finite-mass-qed-absolute-vacuum-window-2026-11-27.json).
