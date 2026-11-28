# Cross the frozen archive histogram and the native decay mesh

The [frozen Hartree check](INDEPENDENT-FINITE-MASS-QED-COLD-HARTREE-SENSITIVITY-BUDGET-2026-11-28.md) traces a narrow cold reaction aperture's mass sensitivity to a 1.67% difference in the archived scalar count. Here test the original massive-daughter 1↔2 production benchmark at x=.05, before adding the separate photon-assisted transfer. Compare each frozen archive on both its own native mesh and the other mesh.

The remap averages occupations with exact overlap lengths in p³, so each species' number and the Pauli caps are preserved. It does not preserve each cell's energy moments, and it changes the representation of a sharply peaked spectrum. For each remapped or native state use the original driver at the same x, recompute the Hartree mass and evaluate the right-hand side with the same radiation correction. Each native replay matches the archived growth clock and mass. Both remaps conserve species number to 3e−15 relative, and the collision source closes the paired number ledgers to 4.2e−17 absolute or better.

| Frozen input | Fine-mesh d log N / d log a | Third-mesh d log N / d log a |
| --- | ---: | ---: |
| Fine | 4.457997 | 2.852987 |
| Third | 5.173667 | 2.189097 |

Holding the archived source fixed before number-preserving remap, the fine/third growth-clock ratio is 1.56 (fine input) or 2.36 (third input). The comparison combines destination discretization with its finite-cell re-representation; it is not a pure interpolation-error estimate. Both effects merit testing before attempting a time-integration tolerance fix. In every arm more than 99.99% of the net scalar source lies below q=1e−5. The Hartree masses on the two destination meshes for a fixed input are close by comparison, so matching the instantaneous mass alone does not explain the clock discrepancy.

Run `OPENBLAS_NUM_THREADS=1 python -W error early_pole_endpoint_stencil_remap_pilot.py` from `kinetics/` and compare the [receipt](results/early-pole-endpoint-stencil-remap-2026-11-28.json). The four instantaneous sources do not isolate the origin of the earlier trajectory difference, guarantee convergence under refinement, evolve the remapped state, or predict a physical photon-assisted rate.
