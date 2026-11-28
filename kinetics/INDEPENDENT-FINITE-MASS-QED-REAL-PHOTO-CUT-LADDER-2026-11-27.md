# Bare real photoabsorption grows as the infrared cut is lowered

Integrate the vacuum reaction $P+\gamma\leftrightarrow b+s$ at
archived scalar momentum 4.1e−6, with frozen charged and scalar
histograms and photon temperature 1. Use a bath-frame **lower**
cut on the real photon energy. This is a different process from
[crossed scalar absorption](INDEPENDENT-FINITE-MASS-QED-REAL-SCALAR-ABSORPTION-2026-11-27.md).
Only one charge assignment is used, with unaveraged spin and
polarization sums. No eikonal or medium counterterm is included
in the numbers in the table.

Use inverse two-body coordinates with fixed incoming daughter
energy and angle to the prescribed scalar mode. The real photon
opens only a narrow distance above the invariant pair threshold.
Resolve that distance with a scaling grid, split the
photon direction near its collinear peak, and split the lab
daughter angle at exact intersections of charged occupation
boundaries and the photon limits. A global unsplit radial
mesh near the opening did not reproduce the low-cut result.
The same chart agrees with the prior restricted comparison
when its photon interval and two charged momentum intervals
are reinstated and its factor of two for conjugate charges is
removed.

This completed ladder uses the fine archive.

| Lab lower photon energy | Unscaled real gain | Unscaled real loss | Signed bare net / Born at g=.005 |
|---:|---:|---:|---:|
| 4e-05 | 132051.2 | 37024.7 | +0.05860619 |
| 0.00016 | 33090.27 | 8944.228 | +0.01489171 |
| 0.00064 | 8317.876 | 2096.925 | +0.003836677 |
| 0.00256 | 2103.724 | 489.7571 | +0.0009953898 |
| 0.01024 | 545.2911 | 91.41976 | +0.0002799183 |
| 0.04096 | 148.9436 | 0.2027809 | +9.173369e-05 |
| 0.5 | 13.53194 | 4.847161e-06 | +8.345618e-06 |

Across the ladder the largest two-mesh relative drift is 8.15e−5.
The restricted-window gain and loss agree with the earlier archive within 1.63e−5.
The 4e−5 products are 1.00123 and 0.98051 times the independent gain and loss soft scales.
The finite parent-support/Planck bound on the remaining unscaled gain is 8.47e−6,
and the corresponding inverse loss tail is zero. Common-temperature net/turnover is
1.05e−16 in magnitude after forming the quantum scalar factors.

To check the normalization independently, multiply the bare
gain and loss integrals by the photon cutoff. The low-cut
limits at fixed photon temperature are $T B/(4\pi^2)$ times
the respective Born gain and loss, with
$B=2(\chi/\tanh\chi-1)$ and
$\cosh\chi=(P\cdot b)/(M_P m_b)$. This follows from the
angularly averaged eikonal radiator and the leading Bose
weight $T/k$. The tested low cutoff is nonzero; its products
need not exactly reproduce these analytic limits.

**Do not add this cutoff-dependent percentage to the Born
source as a physical correction.** Emission, the thermal virtual
term and the finite matching are not included. The large size
and cutoff sensitivity of the bare real contribution make local
soft subtraction preferable to subtracting two separately
rounded large integrals. Even a numerically well-resolved
cut ladder does not determine the sign of the inclusive rate.

Replay from `kinetics/` with the pinned archives:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_hard_photo_complement_reference.py fine
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_hard_photo_cut_audit.py fine
~~~
