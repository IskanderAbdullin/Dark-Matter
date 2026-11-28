# Cross-check the signed elastic number rate in the first light shell

The [tagged elastic calculation](INDEPENDENT-WARM-COLD-LIGHT-BELOW-POLE-ELASTIC-2026-11-26.md)
measured a rest-point rate, not the number moment advanced by
the finite-shell solver. Independently integrate its exact
near-rest (Lphi	o Lphi) collision rate with (p^2dp)
over the whole first shell, imposing (sqrt{s}/Mle.82).
Check directly against the signed number increment of the native
conservative elastic projection at the same archived shape.
Compare the first and fourth matched frames at q16 and q32.

The independent integral changes by about (10^{-5}) relative
or less across two radial and collision-quadrature ladders.
At the first q32 frame it gives (-1.14248	imes10^{-5})
in shell-averaged occupation per Hubble time. The native
((3,2,2,8)) setting underestimates its magnitude by 9–10%
across the two tested q32 frames. At ((9,6,6,16)) the two
methods agree to better than 1% on all four cases, with
the same negative sign in every calculated cell average.
The native operator independently conserves both integrated
species numbers and total dark energy in its represented domain.

This is a local test of a small signed difference between elastic
gain and loss. It neither tests the pole nor sets an error bound
for all momenta, all interaction types, or physical evolution.
Part of the remaining difference can come from native represented
cutoffs and their thermal quadrature conventions.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_cold_elastic_projection_audit_reference.py
~~~

[Diagnostics](results/warm-cold-elastic-projection-audit-2026-11-26.json).
