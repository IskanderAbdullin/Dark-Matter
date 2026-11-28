# Positive and entropy-monotone discrete updates on kinetic-offset shells

Combine the [two-daughter event deposition](INDEPENDENT-TWO-MASS-KINETIC-OFFSET-TRACE-2026-11-27.md)
with a manufactured finite occupation control. Use a shared parent
shell and three different admissible decay directions. Each direction
has its own paired daughter weights on one common momentum chart.
The stoichiometric vector ν removes one parent and adds one scalar
and one fermion in total. It simultaneously adds the separately
supplied parent rest excess to the kinetic ledger.

For count n_i and capacity g_i, use the bin occupation f_i=n_i/g_i.
The Bose and Fermi entropy conjugates are log[(1+f)/f] and
log[(1−f)/f], respectively. Interpolate these conjugates with the
same positive event weights; define A=ν·α. Use a signed quantum
flux from the resulting effective occupancies, rather than
subtracting large nearly balanced products. The flux S has the
sign of A, so the discrete entropy production is S A≥0.

For each event, update the common extent ξ implicitly. The intersection
of all positivity and Pauli inequalities defines a feasible interval.
Along it dA/dξ=Σν_i² s_i″(f_i)/g_i<0. The unique equilibrium point
therefore brackets an entropy-increasing step between the current
state and the equilibrium, regardless of whether the source is
forward or inverse. Capacity inequalities are never imposed by
clipping the result. The fixed-grid energy and number identities
remain linear in the extent.

The fixture uses ordinary and near-rest kinematics. Thermal
equilibrium affinities stay below 3.8e−19 in both. Eight controls
start with perturbed parents, almost-empty daughters, or crowded
daughters close to the fermion ceiling, and take steps up to 50 in
arbitrary event units. Both forward and inverse sequences keep
positive counts and the Pauli ceiling, increase entropy, and
close the split kinetic-rest balance on the small kinetic scale
to below 4.2e−20 in these runs.

This is a manufactured structure-preserving update, not yet a
convergence result for coupled collision kernels. The chart
capacities and event rates are not a microscopic computation.
No physical occupation, tail abundance or cosmological history
is inferred.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_entropy_extent_reference.py
```
[Full receipt](results/two-mass-entropy-extent-2026-11-27.json).
