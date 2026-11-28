# Resolve the active photon spectrum rather than forcing Bose equilibrium

Replace the whole-spectrum Bose fit in the fixed-scale
[photon number budget test](INDEPENDENT-EARLY-POLE-PHOTON-DECAY-COMPARATOR-2026-11-26.md)
with evolving occupations at the hard incident quadrature nodes.
Photons below p/Tref=.5 or above 40 remain unchanged; these spectator
populations contribute their initial number and energy. Fermions
still use the common FD moment target. Every event removes exactly
two photons or two fermions and deposits the same pair energy in
the other sector, with forward and reverse occupation factors.
The Bose plus FD equilibrium test has vanishing net source.

For the hot initial bath, g=.001 and reference exposure .8,
the finer 20-node calculation gives daughter count .12468,
versus .13001 when all photons are forced to a number- and
energy-matched BE shape. The evolving active spectrum differs
from the BE fit to its actual whole-photon moments by about
5.8% in integrated absolute count. The 14-node result gives
daughter .12475 and shape difference 5.9%. The colder bath
has a larger shape difference and a distinct parent-annihilation
initial condition; its finer daughter count is .02923. The
node tests conserve total gauge particle count and energy to
roundoff in the nonexpanding box. Hard phase-space cuts and
thermostat assumptions are not physical scattering boundaries.

This model omits photon momentum exchange, soft processes,
decays and expansion; the BE comparator instead assumes
arbitrarily fast momentum redistribution. The two are
sensitivity cases, not a rigorous bracket on the true bath.

~~~sh
python -W error early_pole_spectral_photon_pair_box_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --order 20
~~~

[Raw photon shape and count/energy ledgers](results/early-pole-spectral-photon-box-2026-11-26.json.gz).
