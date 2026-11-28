# Reflected azimuth stencil on both signs of the native transfer bridge

The existing reduced spinor kernel is invariant under reflection in the
x–z reaction plane. Reflecting the y component conjugates the real-vector
slash representation, with an additional sign for the y polarization,
which drops out of the squared amplitude. An even midpoint azimuth rule
can therefore retain just half the points and average them with equal
weight. This is the same quadrature, not an independent refinement.

Test both variants on six q values per native archive (fresh points just
below and above the energy crossing plus interior points from each band
side), three cosines per q, and two radial/azimuth orders. Both sides have
nonempty supported probes. Across 72 stencil comparisons, the largest
absolute forward/inverse/net difference normalized by gross forward plus
inverse is below 4.0e−14; the support flags agree. This is a local control
on future CPU reduction, not a complete bridge rate, a replacement for
the independent medium mesh, or a physical prediction. The existing
full-stencil calculations were not interrupted or replaced.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error finite_mass_qed_native_transfer_mirror_control.py
```
