# Massive daughter removal and inverse production

The stable massive daughter in [the two-mass calculation](TWO-MASS-BATH-2026-11-25.md) has an unacceptable matter density. Add a phenomenological detailed-balanced interaction with the thermal reservoir:

d f2 / d ln(a) |_bath = (Gamma2/H) <m2/E2>_cell [f_FD(E2/T_bath)-f2].

kappa2=Gamma2/H_ref(M); the actual H includes tracked parent, daughters and bath. The shell average uses the pressure-trace weight. Inverse creation is included, unlike a pure sink. The same reservoir exchanges the corresponding energy, and an independent ODE variable integrates the net fermion number source. The identity N2-Nphi-S=0 and total FRW continuity are checked in every recorded state. The model assumes constant g_bath and vacuum dispersion; it is not a microscopic finite-temperature decay calculation.

## Executed exploratory cases

Windows Python 3.9.7, NumPy 1.26.4, SciPy 1.13.1; single BLAS thread. n=160, m2/M=1/sqrt(2), mphi/M=1e-9, x0=.1, xend=80, parent chemical relaxation/H=100, kappa=.05.

| kappa2 | omega_phi h² | p/T>.1 fraction | final T_bath/T_ref |
|---:|---:|---:|---:|
| .01 | .139446 | .0430794 | 1.025842 |
| .03 | .143090 | .0493784 | 1.019295 |
| .1 | .146817 | .0636423 | 1.015691 |

All three runs passed thermal balance, directional Jacobian, independent integrated number source, Pauli and continuity checks. Final fermion abundances are near the numerical floor (kappa2=.01: N2/Tref³=3.12e-14; .03: -5.82e-18; .1: 9.41e-19). Negative residuals are numerical, NOT negative physical density. The formal rest density assuming all reactions stopped at x=80 for .01 is ~1.03e-5; for the other two cases no precise residual relic density is inferred. The actual massive daughter is unstable in this model and is not counted as constant Delta Neff. The scalar is still overproduced at this width.

[Source](daughter_bath_relaxation.py), [selected actual outputs](results/daughter-bath-2026-11-25.json), [width tuning](tune_daughter_bath.py). Next tune kappa at fixed kappa2, compare independent meshes and include realistic temperature dependence of the parent/daughter bath rates. A physically viable model additionally needs specified interactions and thermal dispersion; these effective rates alone do not prove viability.
