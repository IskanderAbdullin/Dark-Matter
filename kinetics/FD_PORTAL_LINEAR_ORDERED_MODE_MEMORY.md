# Ordered local slow-mode memory on the sampled expanding paths

Interpolate each accepted direct temperature history on log(a/a_start)=0.25–0.55 and independently reconstruct the two-force, local equilibrium relaxation matrix along it. Propagate an initial unit slow eigenvector in a fixed force scale. Compare with the scalar estimate obtained by integrating the instantaneous slow eigenvalue alone. The mode directions change with temperature; the matrix propagation is not assumed to commute. The two phase-space resolutions and coefficient sampling levels agree on the reported components within 1.59e-8.

For the combined history the integrated slow clock is about 1.0035, while the integrated fast clock is about 269.70. The scalar slow survival is 0.36659; the ordered projection on the terminal unit slow direction is 0.37046, about 1.05% higher, with a fast projection of 1.21e-4. The timelike-only case gives similar differences. A very large fast clock therefore does not remove the slow memory in this diagnostic.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_linear_ordered_mode_memory.py \
 eos-range.dat results/fd-portal-direct-expansion-tracking-continuation-2026-11-28.json
```

Use the earlier local EOS normalization. This is an ordered local equilibrium linearization on an interpolated temperature history, not the exact linearization about the nonequilibrium trajectory or a physical abundance result.
