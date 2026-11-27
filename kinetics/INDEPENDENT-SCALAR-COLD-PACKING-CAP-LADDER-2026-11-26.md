# Occupation ladder for a selected linear-power target

The [native near-packed comparison](INDEPENDENT-SCALAR-COLD-PACKING-NATIVE-CLASS-2026-11-26.md)
used two widely separated scalar occupation scales.
Here keep exactly the same half-abundance $86$ eV
scalar setup, smoothed filled-ball shapes, native
transfer calculation and matched CDM controls.
Only rescale the distribution according to
$q_0\propto F^{-1/3}$ to change the packing
occupancy while preserving the particle count.

For an illustrative numerical target, request
$P/P_{\rm CDM}\ge .95$ at six **sampled** modes
from $.05$ to $5\,h\,{\rm Mpc}^{-1}$, at $z=5,20$.
This target is not a measured Lyman-alpha limit.
The smallest ratio for every new arm is at the
highest retained mode.

| Occupation scale | Precision arm | $P/P_{\rm CDM}$ at $k=5$, $z=5$ | Sampled .95 target at both redshifts |
|---:|---|---:|---|
| 100 | compact | .48545 | fails |
| 1,000 | compact | .83869 | fails |
| 5,000 | dense | .94179 | fails |
| 7,000 | dense | .95312 | passes |
| 10,000 | compact / dense | .96192 / .96279 | passes |

The 10,000 compact/dense arms differ by at most
.091% relative in the tested ratio. The crossings
between the sampled 5,000 and 7,000 scales are
diagnostics in this shape family, not an extremal
occupation requirement for arbitrary spectra.
They do not establish power on unsampled modes
or identify any physical production process.
An observational comparison would need the
data, nonlinear effects and likelihood treatment.

~~~sh
python -W error scalar_cold_packing_native_cap_ladder_reference.py "$BASE" "$LADDER" compact 100 1000 10000
python -W error scalar_cold_packing_native_cap_ladder_reference.py "$BASE" "$LADDER" dense 5000 7000 10000
for cfg in "$LADDER"/*.ini; do
    CLASS_EXACTK_POWER_ONLY=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
      "$CLASS" "$cfg" > "$cfg.log"
done
python -W error scalar_cold_packing_native_cap_ladder_audit_reference.py "$BASE" "$LADDER"
~~~

[Exact ratios and selected target checks](results/scalar-cold-packing-cap-ladder-2026-11-26.json);
[compact raw grids](results/scalar-cold-packing-cap-ladder-compact-raw-2026-11-26.json.gz),
[dense raw grids](results/scalar-cold-packing-cap-ladder-dense-raw-2026-11-26.json.gz).
