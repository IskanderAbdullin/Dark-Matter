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
The same exact coordinate maps were tested in the **unmodified**
solver on a common short k range, about 0.09–0.60 h/Mpc.
All three runs have the same present stable density,
$\omega_{\rm stable}=0.11816900374571584$, and do not include
the extra scalar population. At their largest absolute
departure from matched CDM, z=0:

| Native transformed quadrature | Signed power residual | Independently sampled number-dilation residual |
| --- | ---: | ---: |
| 128, S=1 | −0.00951085 | −0.01141092 |
| 128, S=3.81419 | +0.00038124 | +0.00045486 |
| 224, S=3.89953 | −0.00016891 | −0.00020392 |

Results at z=5 are almost identical. At fixed 128 nodes,
relabeling the same unsmoothed input reduces the largest native
departure from 0.95% to 0.038%, without changing the physical
spectrum or density normalization. The ratio of the two signed
columns is about 0.83 empirically in these three cases, consistent
with the earlier causal derivative diagnostic. It is not a
universal correction formula. Numerical pressure differences
persist (present pressure-to-energy ratios differ by about
0.07% between the two 128-node variants); the residual spread
and weak sign changes are reasons not to claim a measured small
physical transfer. No production-era perturbations are evolved.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error class_exact_coordinate_remap_reference.py \
    fine/stable.dat class-configs/stable-short-native128-2026-11-26.ini native \
    --nodes 128 224 320
~~~

The template must match the intended synchronous model and
contain the stable species, density, hierarchy and desired
short k range. The script rewrites only the output/PSD
paths, the coordinate scale and the numerical quadrature.

For the numerical comparison, run CLASS from the working directory
containing fine/stable.dat and the generated native directory.
The two checked-in short configs reproduce the S=1 control
and the matched CDM background:

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/path/to/class class-configs/stable-short-native128-2026-11-26.ini
/path/to/class class-configs/stable-short-cdm-reference-2026-11-26.ini
/path/to/class native/psd_coordinate_remap128_sync.ini
/path/to/class native/psd_coordinate_remap224_sync.ini
~~~

After those four runs:

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error class_exact_coordinate_remap_audit_reference.py \
    native --stable-table fine/stable.dat \
    --cases psd_coordinate_native128_sync psd_coordinate_remap128_sync psd_coordinate_remap224_sync
~~~

The output also records the density and pressure checks. It uses
the existing short native CDM reference cdm_short_altq_sync.
