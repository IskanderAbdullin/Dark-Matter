# Measure scalar spectrum dependence at the unamplified fraction

The [amplified controls](INDEPENDENT-WARM-DIRECT-NATIVE-PILOT-2026-11-26.md)
used scalar fraction 1e−3 to compare distinct shapes.
Use the relaxed-history q64 spectrum and the independent
direct-production spectrum at the **same** smaller
conditional fraction 2.099185e−6. Use the early-start RK
setup and exact native wavenumbers. Center each shape
on its own same-table 1e−10 trace control. This removes
what would otherwise be a potentially misleading
table-dependent numerical offset as the abundance tends
to zero.

The shape contrast is
ln(P_held/P_direct) at the physical fraction minus
ln(P_held/P_direct) in the trace pair. Equivalently,
it is the difference of each shape's own trace-centered
response; both formulas agree numerically.

| k (h/Mpc) | Trace-centered held minus direct, Delta ln P at z=5 | Same at z=20 |
| --- | ---: | ---: |
| .3 | −3.295e−7 | −3.286e−7 |
| 1 | −9.212e−7 | −9.125e−7 |

At k=1 the response corresponds to about 18% of the
direct scalar's own power effect in this comparison.
The k=1 trace table offset is only 3.66e−9. The amplified
standard-tier tangent would suggest roughly −9.54e−7
at z=5; the result above is measured directly at the
physical fraction instead of extrapolated. At smaller
modes, differences fall toward the solver/numerical
floor. Do not treat the table as an observational
sensitivity or a solution for different physical
production abundances, CDM histories or interactions.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_direct_vacuum_native_early_start_reference.py focus early --opacity 2e-7 --templates native
# The direct trace and physical outputs from the earlier exact-mode control are also needed.
for n in heldtrace heldphysical; do
 /path/to/class "early/$n.ini" >"early/$n.log"
done
python -W error warm_direct_vacuum_native_physical_shape_reference.py early \
 results/warm-direct-vacuum-native-pilot-smoothing-2026-11-26.json
~~~

[Two-trace native shape results](results/warm-direct-vacuum-native-physical-shape-2026-11-26.json).

[Native held-table power and startup rows](results/warm-direct-vacuum-native-two-trace-held-raw-2026-11-26.json.gz).
