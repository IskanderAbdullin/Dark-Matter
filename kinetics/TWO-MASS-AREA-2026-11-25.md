# Fast two-mass daughter-energy area — 2026-11-25

Stage completed locally (NumPy 2.3.5, SciPy 1.17.0). This is a geometric kernel, not yet an evolved distribution. With parent mass x and daughter masses m,b, define K=(x²-m²-b²)/2, D=sqrt((K-mb)(K+mb)). At fixed fermion energy r, p=sqrt(r²-b²), scalar bounds are E_±=(Kr±Dp)/b². Rationalize E_-=m²r/(K+D)+D/(r+p). E_- has minimum m at r=K/m, while E_+ is monotone. For b=0 interchange daughter coordinates or use the two-massless case.

Full/empty rectangles are classified in vectorized form. For mixed cells split at intersections with each scalar energy boundary, integrate the appropriate clipped height using analytic primitives. For r=b cosh t, integral dr/(r+p)=t/2+exp(-2t)/4, integral p dr=b²(sinh(2t)-2t)/4. Difference formulas use log1p/expm1 and a near-rest rapidity series. Extended precision on mixed cells minimizes cancellation. No approximate thermal dispersion is included.

**Actual local checks:** 8509 nonzero-width rectangles, including targeted kinematic boundaries and near threshold, compared to independent adaptive quadrature; maximum absolute area-fraction discrepancy 4.38e-9, daughter exchange 2.06e-11, 600 subdivision checks 2.17e-16. Adaptive quadrature emitted 121 warnings in extremely thin/near-threshold cells; comparison is not a certification of their reference precision. The random test has broad but finite coverage. Thirty 480×480 matrices took 0.297 s on this local machine (not a universal performance promise). Threshold and nonnegativity checked. Inputs whose energy edges become identical at float64 precision give zero measure; this requires explicit momentum-grid error accounting in the evolution solver. No post hoc spectrum adjustment.

- [Implementation](two_mass_area.py)
- [Reproduction](test_two_mass_area.py): from kinetics, run OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python test_two_mass_area.py.
- [Raw results](results/two-mass-area-2026-11-25.json)
- Previous independent reference: [two-mass geometry](TWO-MASS-GEOMETRY-2026-11-25.md).

Next: integrate the matrix in the energy-closed solver with both massive shell energies and traces; retain number conservation and distinguish massive daughter relic abundance from relativistic Delta Neff. Explore actual production couplings and thermal mass/width consistency; do not infer physical viability from this kernel.
