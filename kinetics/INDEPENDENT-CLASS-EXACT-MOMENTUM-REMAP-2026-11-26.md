# Exact momentum reparameterization for native numerical sampling

The [massive moment check](INDEPENDENT-CLASS-MASSIVE-MOMENT-DILATION-2026-11-26.md)
shows that simply increasing the number of raw momentum nodes
does not suppress all derivative errors. Select an alternative
grid without smoothing the archived distribution or replacing
the CLASS derivative. With a dimensionless positive scale $S$,
relabel the momentum by $q'=q/S$, tabulate $f'(q')=f(Sq')$,
and set $T_{\rm ncdm}$ to $S$ (relative to the original 1).
The same physical momentum maps to the same occupancy.
The effective mass parameter in the integrand changes to
$(am/T_{\gamma,0})/S$, so all energy-weighted integrals change
by $S^{-4}$ in dimensionless coordinates. The external factor
$T_{\rm ncdm}^{4}$ restores the physical density and stresses
exactly. The corresponding number-weighted factors are
$S^{-3}$ in the integral and $S^3$ outside. This is a coordinate
change, not a warmer or cooler physical distribution.

Apply the unmodified indefinite trapezoidal quadrature
$q'=1/t-1$ after this relabeling. Scan the scale using only
external background energy, pressure and their dilation
identities across six dimensionless mass-to-temperature
ratios. No power or perturbation outputs are used to choose $S$.

| Nodes | Worst screened residual, S=1 | Selected S | Worst screened residual after remap |
| ---: | ---: | ---: | ---: |
| 128 | 0.01141 | 3.81419 | 0.0009454 |
| 224 | 0.00628 | 3.89953 | 0.0002664 |
| 320 | 0.004265 | 2.97697 | 0.0000923 |

The continuum covariance tests for energy, pressure and
their derivative moments agree to better than 4e-15.
The selected points minimize a discrete scan, not a
certified quadrature tolerance; untested epochs or the
compiled solver's finite derivative may behave differently.
This artifact prepares reproducible short-k native configs.
It does not give a native transfer result or evidence about
the actual abundance, observations or convergence of CLASS.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error class_exact_coordinate_remap_reference.py \
    fine/stable.dat stable-ncdm-template.ini native \
    --nodes 128 224 320
~~~

The template must match the intended synchronous model and
contain the stable species, density, hierarchy and desired
short k range. The script rewrites only the output/PSD
paths, the coordinate scale and the numerical quadrature.
