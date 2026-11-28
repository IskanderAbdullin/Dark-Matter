# The cold reaction aperture amplifies the archived Hartree number discrepancy

The [pointwise intervention test](INDEPENDENT-FINITE-MASS-QED-COLD-SENSITIVITY-2026-11-28.md) identified the frozen scalar pole mass as the largest contributor to the difference in raw cold affine coefficients. This check reconstructs where that mass difference comes from in the two archives. It does not evolve either archive with a common number density.

Use the scalar shell counts and the existing repulsive quartic gap equation,
```text
m² − m0² = (lambda/2) integral dN(q)/sqrt(q²+m²).
```
The saved count of every shell equals its constant occupation times the exact spherical volume. With the archived momentum nodes, the archived mass solves the same gap equation to 3e−13 relative or better. Replacing each node by direct quadrature across its own shell changes the mass by only 2.85e−7 (fine) and 1.26e−7 (third) relative. Two independent shell quadrature orders agree within 2e−13 in the reconstructed mass.

The scalar number differs by 1.6748% between the two snapshots; their exact-step Hartree masses differ by 0.55513%. More than 0.9999982 of either archived tadpole comes from q<1e−5. Independently normalizing the fine scalar histogram to the third's total number while retaining its shape reduces the mass difference to 3.86e−7 relative. The actual pole difference is thus a population-normalization effect in this frozen comparison, not an artifact of approximating the tadpole at shell nodes. The numerical elasticity d log(m)/d log(N) is about 1/3 on each archive, confirmed by an independent density variation.

The saved time slices show where to test remeshing first: the yield difference is only about 0.084% at x=.02, rises to 2.52% at x=.04, and is 1.67% at x=.05. Sparse checkpoints establish a useful interval to inspect, not the time or numerical mechanism of the divergence.

Holding the charged distributions and photon dispersion fixed, measure d log(A)/d log(m) at q=3.5e−6 and 8.2e−6 with two symmetric mass increments. The local raw forward sensitivities are 274.8 and 250.0. Combined with the Hartree elasticity, a hypothetical 1% *local raw forward* tolerance requires the count to be controlled near 1.09e−4 and 1.19e−4, respectively, along this one perturbation. This is a numerical target on the frozen Hartree contribution. It does not include the explicit scalar occupation factors in the full collision source, variations of the photon pole, or a self-consistent evolution across archives.

Replay `OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_cold_hartree_sensitivity_budget.py` from `kinetics/` and compare the [receipt](results/finite-mass-qed-cold-hartree-sensitivity-budget-2026-11-28.json).
