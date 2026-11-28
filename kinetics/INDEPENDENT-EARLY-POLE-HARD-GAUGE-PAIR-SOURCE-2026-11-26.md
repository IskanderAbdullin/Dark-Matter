# A populated shared photon bath opens a hard fermion pair source

Consider just the central part of the tree-level
[Breit–Wheeler channel](https://doi.org/10.1103/PhysRev.46.1087),
γγ↔χ_dχ̄_d, in the shared Abelian bath. Use two transverse
massless photon states with thermal Bose distributions, each
incident energy above 0.5 times the bath temperature, s above
its square, and |cos θ*|≤0.8 in the pair center-of-mass frame.
The final cut avoids the collinear angular part of the vacuum
amplitude. Integrate the central cross section analytically in
angle and check it against direct quadrature and the total
vacuum formula.

For the common-temperature FD daughter target with unresolved
chemical deficit, forward minus inverse is forward times
1−exp(2μ_d/Tbath) by detailed balance. Bound outgoing blocking
within this vacuum cut using the largest allowed FD occupation.
The identical incident factor 1/2 and two particles produced
per event cancel in the count rate.

At g=.001 the restricted initial net source in the retained
ratio-.06 benchmark is about 0.272 Tref³ per reference log
expansion step at x=.05, and 0.381–0.383 at x=.10. The snapshot
number gaps to the fixed-parent inversion gate are about 0.33
and 0.16 in Tref³ units. The source varies as roughly g⁴ over
the examined weak-coupling interval. Comparing an instantaneous
source to this gap is not an equilibration time: the photon
population, fermion shape, inverse rate, parent reservoir,
and temperature may evolve. This vacuum hard subset is not a
rigorous lower bound on real in-medium chemical equilibration.

~~~sh
python -W error early_pole_fermion_population_gate_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz > population-gate.json
python -W error early_pole_hard_gauge_pair_source_reference.py population-gate.json
~~~

[Raw hard subset](results/early-pole-hard-gauge-pair-source-2026-11-26.json).
