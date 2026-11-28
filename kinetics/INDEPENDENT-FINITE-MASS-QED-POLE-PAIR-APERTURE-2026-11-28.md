# The frozen timelike collective poles cannot open a free pair cut

The [transverse plasmon](INDEPENDENT-FINITE-MASS-QED-FROZEN-TRANSVERSE-PLASMON-GATE-2026-11-27.md) and [longitudinal Coulomb pole](INDEPENDENT-FINITE-MASS-QED-LONGITUDINAL-SUM-2026-11-28.md) are above the phase-velocity-one line in their tested timelike domains. The statistical Landau cuts computed from a free charged line scattering into another identical-mass on-shell line have spacelike support. The other possible free one-loop cut creates or annihilates a same-species charged pair; its invariant threshold is K² ≥ 4m², regardless of Pauli occupancies. Comparing just ω with 2m at k>0 would be the wrong threshold.

In the frozen hard-loop dispersion, the entire candidate transverse and timelike longitudinal branches are below that threshold. An analytic bound follows without an extrapolated rate. Write β=k/ω, v=p/E, E=√(p²+m²); on these timelike branches 0≤β<1.

**Transverse.** Its invariant is Π_T(β). The radial response can be integrated by parts as an integral of f times the p derivative of p³ h(βv)/E. The angular h has a power series in its argument squared with all coefficients positive, and both p³/E and v increase with p. The response is therefore nondecreasing in β for any nonnegative isotropic occupation with finite moments. At β=1 the p primitive equals (pE−m²atanh v)/2, whose derivative is p²/E. Integration against f yields the exact upper bound within this model. It is the asymptotic transverse invariant, not a prediction for the large-momentum interacting theory.

**Longitudinal.** Expand the real density response in β², Π_00 = Σ_(n≥1) A_n β^(2n). The radial integration-by-parts weight for A_n is proportional to

```text
E v^(2n) [1 − (2n−1)/(2n+1) v²] f(p),
```

which is positive and decreases pointwise with n. The longitudinal root has k²=Π_00; its invariant is (1−β²) Π_00/β². The coefficient of β^(2n), n≥1, in that product is A_(n+1)−A_n≤0. Hence the invariant never exceeds A₁, the square of the plasma gap at k=0. This argument does not need a numerical decision close to the light cone.

| Frozen endpoint | maximum longitudinal invariant mass | maximum transverse invariant mass | lowest free pair threshold |
| --- | ---: | ---: | ---: |
| fine | .001663861 | .002035880 | .006000000 |
| third | .001665235 | .002037563 | .006000000 |

The upper masses are below 34% of the minimum exact free pair threshold. Independent positive radial integrals confirm the jump-form bounds to 6.67e−16 relative; a 96-point β sweep per spectrum and the positive coefficient ladder provide numerical regression guards. The kinematics considered are the free one-loop cut structure and a frozen dispersive reference; a self-consistent interacting pole near a changed threshold was not computed. For the decomposition of free one-loop vector cuts see [DeRocco et al., Appendix B](https://arxiv.org/abs/2201.05167).

Run `python -W error finite_mass_qed_collective_pole_pair_threshold_reference.py` from `kinetics/` and compare the [receipt](results/finite-mass-qed-collective-pole-pair-threshold-2026-11-28.json). This is a specific reason that neither the free scattering noise nor a free-pair gain/loss cut can set the occupation of these timelike poles. Widths or occupation changes due to other interactions, and an independent collective initial condition, remain outside this check.
