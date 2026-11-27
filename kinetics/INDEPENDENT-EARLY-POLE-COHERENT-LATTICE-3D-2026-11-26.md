# A small three-dimensional coherent-field control

Extend the earlier [2D reference](INDEPENDENT-EARLY-POLE-COHERENT-LATTICE-2D-2026-11-26.md)
to a periodic cube with three resonant Cartesian waves of fixed
physical dimensionless wavenumber and matched relative seed
weights and phases. The quartic real field is evolved by
spectral derivatives, a de-aliased force and symplectic steps.
A smaller-wavenumber matched input acts as a negative control.

Across 32³ and 48³ and time steps 0.06 and 0.045, the early
amplitude slope is 0.03610–0.03624 per homogeneous oscillator
time, near the coherent analytic reference. At the tΩ=360
readout, runs on 32³ and 48³ with the same time step give
sampled maximum inhomogeneous energy fractions near 0.330,
and quadratic resonant-band power near 0.287–0.289 of the
total energy. The maximum relative energy drift is below
7.4 × 10⁻⁴. The stable seed control stays below 3 × 10⁻¹¹
in the inhomogeneous energy diagnostic.

This tests an aligned and artificial initial state, over a
short interval in a small box. It does not resolve the broad
separation between the actual cold momentum peak and the
coherent resonant band, or validate an early expanding dark
spectrum. Inhomogeneous and quadratic energies are separate
diagnostics, not particle-number observables in this model.

~~~sh
python -W error early_pole_coherent_quartic_lattice_3d_reference.py
python -W error early_pole_coherent_quartic_lattice_3d_reference.py --cells 48 --step .045 --end 360
python -W error early_pole_coherent_quartic_lattice_3d_reference.py --stable --end 360
~~~

[Raw 3D traces and controls](results/early-pole-coherent-lattice-3d-controls-2026-11-26.json.gz).
