# Complete the vacuum radiative decay, including all hard photons

The [previous absolute window](INDEPENDENT-FINITE-MASS-QED-ABSOLUTE-VACUUM-WINDOW-2026-11-27.md)
used on-shell massive charged fermions, the MS-bar Yukawa coupling
at the parent mass, and virtual plus unresolved real photons. Here
we also integrate all the resolved real photons, up to the exact
parent-rest endpoint omega/M = [1−(m+ms)²/M²]/2. No final-state
occupation or other thermal weight is applied.

An independent invariant chart integrates over the daughter-photon
invariant t = (b+k)²/M² first. At fixed t the photon energy fraction
r = omega/M ranges between k* (E*−p*) and k* (E*+p*) in units M=1,
where k*=(t−m²)/(2√t), E*=(1+t−ms²)/(2√t) and p* is the parent
momentum in the daughter-photon center of mass. Intersect this
interval with the chosen hard band and split at both band-boundary
intersections. The squared amplitude uses the independent invariants
P·k = r, b·k = (t−m²)/2, and P·b = (1+t−ms²)/2−r. This chart and the
direct photon-energy/polar-angle chart agree on the integrated hard
real coefficients to less than 7.5e−14 absolute in these benchmarks.

On both archived mass choices, vary the artificial soft/hard
separator from .001 to about .353 in omega/M. At .01 the unresolved
plus virtual coefficient is approximately −11.168 and the hard real
coefficient is about +10.748. Their sum does not change with the
separator (maximum spread 1.1e−13). Integrating the matched recoil
directly to the endpoint gives −0.41960951 (fine) and −0.41960906
(third), with absolute agreement better than 2.4e−10 to the separated
computation. The finite photon-mass replay at lambda/M=1e−8 remains
about 3.2e−7 from its regulator-free result.

All coefficients above multiply g²/(4 pi²) in the rate relative to
Born; g=.005 would give a complete vacuum rate correction of about
−2.66e−7. The much larger negative correction from the previously
reported cap alone is not the inclusive vacuum correction.

A separate infrared-finite first-moment test uses (t−m²)/2 for the
loss of scalar energy relative to Born, and r for the emitted photon
energy. In units M times g²/(4 pi²) per Born decay, the mean energy
changes are +0.97491 (photon), −0.13880 (scalar), and −0.83611
(charged daughter, by energy conservation). The photon and scalar
finite-band energy moments both match the reversed invariant
integration. The virtual channel supplies no recoil shift relative
to the Born energies. These are rest-frame vacuum averages, not a
medium source, eventwise energy distribution, or cosmological rate.

Replay from kinetics with the two spectrum archives:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_inclusive_vacuum_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-inclusive-vacuum-2026-11-27.json).
