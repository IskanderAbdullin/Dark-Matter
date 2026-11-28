# Different fermion populations can shift the decay aperture

Equal charges alone do not make the two fermion asymptotic self-energies
equal when their distributions differ. For the leading lightlike
massless-loop proxy, each species contributes
g²/π² times the integral dp p of its own symmetric FD occupation
plus the common photon occupation. This follows from the QED
weights in Table 4 of [Czajka and Mrówczyński](https://arxiv.org/html/1410.5117),
with the hard lightlike mass equal to twice the soft mass parameter.
It reproduces g²T²/4 in the common equilibrium control. This is a
kinematic proxy for hard fermion poles with small bare masses; the
full finite-mass hard self-energies have not been matched.

For the archived x=.05 shapes, the parent-minus-daughter squared-mass
difference from this correction is about .068 g² T_ref² on both grids.
The shared photon term cancels from that difference and changes the
aperture width slightly, not its energy center. At g=.005, the
species-dependent part shifts the parent-energy center by .1447 of
its unshifted half-span on both grids. A common shift leaves the center
unchanged. At g=.01 the same frozen mass prescription shifts it
by .579 half-spans.

To distinguish the shift from interpolation, split the aperture at
both fermion momentum grids and integrate the seed f_P(1−f_b) and the
stimulated difference f_P−f_b. Repeat with constant cells and
positivity- and number-conserving minmod reconstruction in u=p³/3.
Also cross the two occupation grids with the two scalar Hartree poles
instead of comparing only mismatched default poles. In the represented
tests at g=.005, the normalized stimulated aperture increases by about
.07% to 17%, depending on grid and pole. With each grid's own pole
and constant occupations, the increases are 4.2% and 13.6%.
Thus the center displacement is robust within the assumed dispersion,
but the occupied aperture response is not numerically converged.
The sign stays positive on this particular frozen x=.05 snapshot;
there is no inference about the later sign crossing.

No physical width, matrix-element correction, scalar Hartree feedback,
cosmological trajectory, or chemical threshold is calculated here.
Reproduce from kinetics/ with
early_pole_frozen_flavor_asymptotic_mass_reference.py.
Raw crossed matrix: results/early-pole-frozen-flavor-asymptotic-mass-2026-11-26.json.gz.
