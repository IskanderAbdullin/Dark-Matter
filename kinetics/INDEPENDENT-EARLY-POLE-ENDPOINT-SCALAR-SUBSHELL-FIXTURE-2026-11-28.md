# Isolate the cold scalar shell grid from the charge grid

The two-subcell continuation changes the parent and fermion resolution but keeps the scalar grid unchanged. For a separate control hold both charged charts fixed at two subcells. Geometrically subdivide only existing scalar shells fully at or below q=1e−5, preserving all native edges and keeping the native geometric-center convention in each new shell. For the frozen control, copy each original scalar occupation into all its children; the q²dq weights then preserve the total and each original shell's number. This control does not preserve exact scalar energy across the remap, and the refined Hartree mass is recomputed from the actual new nodes.

On the native x=.05 cards, the one-subdivision wrapper reproduces the old chart's short evolution and final state exactly. The charge nodes and weights are byte-identical across all scalar resolutions. A separately constructed thermal occupation passes the microscopic detailed-balance test; exact event energy projections remain below 5e−16 relative. At the frozen occupations, the instantaneous *net* rate normalized by the current scalar number, including the prescribed expansion-clock factor, is:

| Native grid | Scalar shells, unrefined / doubled / quadrupled | Local log-clock number slope at these three resolutions |
| --- | --- | --- |
| Fine | 67 / 93 / 145 | 4.41827 / 4.40122 / 4.39696 |
| Third | 97 / 134 / 208 | 3.67048 / 3.67501 / 3.67661 |

The fine slope drops while the third slope rises. The sensitivity in this frozen control is small but nonzero and cannot establish the sign or magnitude of a change in the full evolving spectrum. The result also does not constrain the accuracy of one particular cold photon aperture. The opt-in startup driver supports independent empty-daughter trajectories on each scalar chart; those trajectories are not included in this result.

From `kinetics/`, reproduce the [receipt](results/early-pole-endpoint-scalar-subshell-fixture-2026-11-28.json):
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error early_pole_endpoint_scalar_subshell_fixture.py
```
