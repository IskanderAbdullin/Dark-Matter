# Use reflection to keep exactly the same even-order azimuth quadrature

In the [isotropic timelike-pole radial integrator](finite_mass_qed_timelike_pole_batched_absorption_reference.py), the photon and scalar momenta lie in the x–z plane. Pair each midpoint azimuth φ with 2π−φ, which reverses only the y components of the two charged momenta. In this code's Dirac basis, γ² and σ² are purely imaginary and the other matrices used in these spinors are real. Reflecting y sends the spinors and slashed real charged momenta to their complex conjugates. The spin-summed matrix with a real polarization other than y conjugates; the y-polarized matrix additionally acquires an overall sign. Its squared norm is unchanged in either case. All energy denominators and isotropic statistical factors are also unchanged, including the optional continuum occupations.

With even full order N, the average can therefore evaluate only the first N/2 azimuths and average them with weight 2/N, rather than averaging all N at 1/N. This is **not** using a lower-order quadrature or imposing a new angular cut. It applies to the existing real linear transverse and longitudinal basis and isotropic unpolarized source. It does not justify the same optimization for general complex coherent polarization or anisotropic occupations.

The [fixture](results/finite-mass-qed-timelike-pole-mirror-fixture-2026-11-28.json) repeats both calculations for each of the two native archives. It checks twelve nonempty radial cases spanning cold, past and future transverse transfer and the longitudinal control, two complete cosine integrals, and an evolved scalar-refined state with smooth charged occupations. The radial coefficients agree below 2e−14 relative to the largest coefficient. The measured isolated radial calls generally use about half the CPU time; this is a local sample, not an end-to-end speed claim. Full-mesh continuation remains the default. The opt-in cold-band driver records its `--mirror-azimuth` choice and rejects mixing it with the legacy raw files in a resumed run.

From `kinetics/`:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_timelike_pole_mirror_fixture.py
```
For a new cold-band invocation use the ordinary full even-order mesh and add `--mirror-azimuth`, with a clean output path. This change leaves the physical kernels, scalar panels, radial quadrature, and the acceptance threshold unmodified.
