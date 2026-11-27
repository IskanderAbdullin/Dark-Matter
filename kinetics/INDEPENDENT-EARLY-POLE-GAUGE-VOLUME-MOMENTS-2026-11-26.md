# An energy-only radial projection can miss cold-event momentum

The [gradient-energy aperture](INDEPENDENT-EARLY-POLE-CONSERVING-GAUGE-EVENTS-2026-11-26.md)
is exactly momentum conserving at the continuous event. Keeping only
a positive linear interpolation of energy onto neighboring radial
nodes does not preserve the on-shell vector momentum of that event.
For the hard convex dispersion, the masses differ between parent and
daughter; their independent energy interpolations change their
conditional mean momenta by different amounts. The resulting absolute
error is tiny on the hard scale but not necessarily on the soft scalar scale.

We use representative linear or volume centers in the two archived
momentum meshes. These are not a measured error of the original solver.
For the frozen g=.005 aperture at q/T=4.1e−6, the energy interpolation
defect is about 2.2e−16 T while the average vector defect is 2.8–5.1%
of the scalar momentum, depending on grid/centers. For roughly 1.6–3.3%
of the geometric parent-aperture length, even readjusting the angle
with the two interpolated mean magnitudes fails the momentum-triangle
test. At q/T=1e−6 the mean vector defect rises to about 12–21% of q.
These are geometric fractions and neither predicted count losses
nor a statement about total vector momentum in an isotropic volume.

**Separate finite-volume moments.** A different representation keeps
the subcell event and accumulates its signed count, energy and two
independent vector components into the physical cells separately.
The support is split where either fermion changes cell. A sine change
of variable removes the integrable angular endpoint behavior.
Across both grids and three scalar momenta, the frozen summed transfer
conserves the respective counts, energy and vector momentum to 2e−15
per unit absolute flux. The independently sampled quantum balance
at common temperature and equal fermion chemical potentials closes
to 1e−15 relative. These separate moment deposits do not construct
an evolving admissible occupation: a subsequent dynamic closure
must reconstruct Pauli-allowed positive shapes from the moments.

Both results intentionally omit the in-medium vertex and use a shared
positive geometry weight. The recovered vector legs form a reference
for a future mesh, not a production spectrum or cosmological
prediction. Reproduce from kinetics/ with
early_pole_gauge_projection_momentum_reference.py and
early_pole_gauge_volume_moment_decay_reference.py.
Raw output: results/early-pole-gauge-volume-moments-2026-11-26.json.gz.
