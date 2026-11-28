# Local adiabatic forcing projects strongly onto the slow FD mode

At common temperature and zero dark chemical parameter, compute the independently adiabatic changes of the two entropy forces per log expansion. Check the derivatives by finite differences of the separate entropy-conserving drift. Project that forcing onto the finite tabulated bath's local timelike and elastic-plus-timelike eigenmodes, then solve the algebraic balance with frozen coefficients. This diagnostic is not a time-integrated lag prediction.

At T=0.55 GeV, the combined fast rate is about 190 times the expansion rate, but the slow rate is only 0.98 times it. The algebraic chemical displacement is -0.202, almost all in the slow mode. That is too large to infer an actual displacement by linearization. At 0.68825 GeV the fast and slow rates are about 1928 and 7.25 times expansion, and the algebraic chemical displacement is -0.0194, dominated by the slow response -0.01947 rather than the fast contribution +7.35e-5. A large fast rate alone does not establish that the bath tracks equilibrium in the chemical direction.

The direct source and projected rates are checked at two phase-space orders over eleven temperatures. The largest scaled adiabatic forcing derivative discrepancy is 3.43e-10, and the largest relative cross-mesh lag discrepancy is 4.42e-6. Large algebraic responses are explicitly labeled as unsuitable for the small-force diagnostic, not predictions of a physical frozen state.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_expansion_forced_modes.py eos-range.dat
```

Use the same authors' tabulated EOS normalization as in [the closed bath benchmark](FD_PORTAL_TABULATED_BATH_ENTROPY_AND_MODES.md). No final abundance or decoupling threshold is inferred.
