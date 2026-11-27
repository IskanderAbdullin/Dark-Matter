# Strict hard interspecies elastic heat

The earlier [projected entropy mixing](INDEPENDENT-EARLY-POLE-ARCHIVED-COLD-ENTROPY-MIXING-2026-11-26.md)
can redistribute energy by construction but does not determine its
microscopic direction. Here use a specific gauge process, elastic
`P+b ↔ P+b`, and retain only the reversible vacuum-hard subset with
incoming `s/T_ref² ≥ 1`, transfer `τ ≥ .15`, and all four momenta
inside a finite radial box. The companion transfer band `[.10, .15]`
is recorded separately. In the massless vacuum limit `τ = |t|/s`.
The [distinct-fermion Born formula](https://cds.cern.ch/record/2712222/files/scoap3-fulltext.pdf)
(Equations 2 and 4) gives `dσ/dcosθ = g⁴[1+(1-τ)²]/(16π s τ²)`
in that limit. There is no screening inside this intentionally hard box.
Initial and final fermions are distinguishable by flavor, and the
Dirac populations have four particle/antiparticle/helicity states
per flavor in the isotropic measure. The factor one half is instead
for integrating forward minus reverse over the full reciprocal
domain. This avoids counting each directed energy transfer twice.

## Frozen snapshots

Both outgoing archive edges and hard boundaries are resolved in
azimuth; for piecewise constant occupations the azimuthal first
energy moment is analytic. The incident grid resolves the largest
occupation changes around the production aperture. Each event uses
the full Pauli-weighted forward and reverse flux and deposits the
opposite energy in the two populations. The uniform-temperature
witness with separate chemical potentials has integrated daughter
heat consistent with roundoff, below
`2e-14` in the reference normalization. The hotter parent witness
heats the daughter in every transfer band.

At `g=.005`, the [.45, 8] radial box gives about 30 in both archives,
in `T_ref⁴` per reference log-scale interval, for `τ ≥ .15`.
Four quadratures (including a low-order split at every incident
occupation edge) give 30.34–30.84 on fine and 29.81–30.82 on third.
Enlarging the box to [.30, 10] gives 29.80 and 30.20 at coarse
resolution. The band [.10, .15] adds about 9.5–9.8 in the [.45, 8]
box. Do not extrapolate this band sum to zero transfer. In the
same radial box as the live test, [.45, 5], the frozen result is
23.78 and 23.99. In contrast, the previous imposed entropy mixing
initially sent energy from the daughter to the parent in its own
active box: the toy mobility is not a proxy for the microscopic sign.

An exact vacuum finite-mass repeat uses the same CM angular cut
`τ = |t|/(4 p_cm²)`, so that the angular cut, not the exact
invariant t cut, is fixed. With `m/T_ref = .05, .003` this changes
the coarse signed heat from 30.4666 to 30.4544 (fine) and from
30.8164 to 30.8037 (third), about 0.04% in both cases. The generalized
formula with zero masses reproduces the original massless integral
to numerical precision. The small mass correction is below the
variation among the tested quadratures and says nothing about
collective medium effects.

## Conservative live projection

The packet quadrature splits the azimuth across both outgoing
radial grids, then uses the same species-dependent energy gradients
as the cold scalar and photon-pair solver. A final daughter radial
coordinate is remapped so that incoming pair variational energy
equals outgoing pair variational energy, and the outgoing momentum
triangle is checked. The same packet represents gain and loss,
with separate flavor counts and total variational energy preserved.
Both entropy production and equal-temperature detailed balance
are independently checked. The dressed remap is a numerical
projection of massless vacuum weights, not a new finite-density
cross section. Neither a continuum soft collision operator, Compton
scattering, same-flavor channels, nor flavor-changing annihilation
is included. No final relic or CLASS input follows.

The live projection initially sends 23.11 (fine) and 23.56 (third,
refined static) to the daughter; the common-temperature moment
residual was below 2e-12. On the fine model at reference exposure
`8e-5`, adding this elastic subset changes the scalar gain from
`6.2342e-6` to `6.4086e-6` without photon pairs (+2.80%) and
from `6.0151e-6` to `6.1918e-6` with photon pairs (+2.94%).
Both use the same initial state and live window; these changes
are not a measurement of full transport.

Run from `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python early_pole_hard_distinct_fermion_elastic_heat_reference.py --level 0
python early_pole_hard_distinct_fermion_elastic_heat_reference.py --level 0 --massive
python early_pole_hard_distinct_fermion_elastic_heat_reference.py --level 1 --without-fixtures
python early_pole_hard_distinct_fermion_elastic_heat_reference.py --level 0 --bounds .3 10 --without-fixtures
python early_pole_archived_cold_elastic_packets_reference.py --grid fine --exposure 8e-5
python early_pole_archived_cold_elastic_packets_reference.py --grid fine --exposure 8e-5 --strength 0
```

Raw: `results/early-pole-hard-distinct-fermion-elastic-2026-11-26.json.gz`.
