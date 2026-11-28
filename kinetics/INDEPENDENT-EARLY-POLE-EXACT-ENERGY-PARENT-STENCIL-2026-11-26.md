# Energy-conserving live-parent stencil

The [resolved parent calculation](INDEPENDENT-EARLY-POLE-RESOLVED-PARENT-2026-11-26.md)
uses exact geometric overlaps and conserves counts, but cannot assume
that its representative shell energies add up in each reaction. This
independent, positive stencil makes that equality exact in the
discrete collision operator. It retains a massless fermion daughter
and an evolving common positive scalar Hartree pole.

For each scalar node and fermion daughter cell, integrate the exact
geometric width across the daughter's aperture. At the two neighboring
parent shell-energy nodes choose barycentric weights that reconstruct
the sum of scalar energy and daughter shell-mean energy. In each
positive forward event deposit one scalar and one daughter and remove
the paired parent weights; reverse all increments for the reverse
event. Forward and reverse parent factors are weighted geometric
means of the occupied and vacant parent fractions respectively.
This makes the Pauli bounds inward at each occupied receiving node.
At a common temperature, the log occupied/vacant ratio of each
parent node is affine in energy; the weights therefore preserve
discrete detailed balance exactly. Scalar and fermion daughters use
their usual Bose and Fermi factors at the corresponding energies.

The interpolation is a discrete closure: there is no claim that
geometric interpolation of arbitrary parent occupations represents
the actual within-shell microphysics. Channels outside the spanned
parent nodal energies are omitted. The omitted fraction of the
**thermal seed** is checked along each trajectory, but this is not a
bound for arbitrary nonequilibrium high-energy parents. The parent
starts FD at the nodes rather than cell-averaged FD, another small
difference from the original resolved calculation.

The static collision checks cover discrete energy reconstruction,
detailed balance, instantaneous total energy at a generic state,
the Pauli face of a parent cell, and a fixed-gap directional derivative.
Across the first two trajectories the normalized instantaneous energy
residual is at most 5.5 × 10⁻¹⁶, and the tested thermal seed omission is
below 10⁻⁵. At x=0.05, coarse/refined scalar yield is 0.000367867 and
0.000350663. On the refined grid this is 1.57% below the matching
cell-overlap method (0.000356239). A third grid gives 0.000350976, a change of +0.089%,
with omitted thermal seed at most 2.37 × 10⁻⁵. This is still a finite-time
early reference, not a freeze-out or observed abundance calculation.

~~~sh
python -W error early_pole_exact_energy_parent_stencil_reference.py --static-only
python -W error early_pole_exact_energy_parent_stencil_reference.py
python -W error early_pole_exact_energy_parent_stencil_reference.py --scalar 65 --fermion 90 --parent 72
python -W error early_pole_exact_energy_parent_stencil_reference.py --scalar 95 --fermion 125 --parent 108
~~~

Raw occupations: [coarse](results/early-pole-exact-energy-parent-stencil-coarse-2026-11-26.json.gz), [refined](results/early-pole-exact-energy-parent-stencil-refined-2026-11-26.json.gz), [third grid](results/early-pole-exact-energy-parent-stencil-third-2026-11-26.json.gz).
