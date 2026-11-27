# Native linear power for two counterfactual near-packed scalar spectra

This is a fixed-abundance check downstream of the
[cold phase-space capacity](INDEPENDENT-SCALAR-COLD-OCCUPATION-PHASE-CAPACITY-2026-11-26.md).
It does not prove that any microphysical source
generates these shapes or that minimizing the
ballistic path also maximizes linear power.

Put 86 eV scalars at $\omega_\phi=.060$ and
an independently imposed CDM component at $\omega_c=.060$;
compare with pure CDM at $\omega_c=.120$.
Keep the remaining parameters fixed. For the scalar
input smooth the ideal occupied ball with a radial
logistic edge $f(u)=1/[1+e^{(u-1)/\epsilon}]$,
where $u=q/q_0$. Restore the number by rescaling
$q_0$ using the integral of $u^2 f(u)$.
Representative geometric caps are 1 and 200,000.
The two combined precision tests change the edge
width (.045 to .0225), hierarchy (28 to 44),
native momentum nodes (96 to 192) and perturbation
tolerance. They do not isolate the source of drift.

Native CLASS 3.3.4, synchronous RK, no ncdm fluid
approximation, exact stored modes, and separately
matched CDM at each tolerance arm. The exact
power-only output patch available in this PR can
suppress unwanted perturbation tables; it does not
change the integration.

~~~sh
python -W error scalar_cold_packing_native_class_reference.py \
 results/scalar-cold-occupation-phase-capacity-2026-11-26.json "$NATIVE_OUT"
for cfg in "$NATIVE_OUT"/{compact,dense}-*.ini; do
    CLASS_EXACTK_POWER_ONLY=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        "$CLASS" "$cfg" > "$cfg.log"
done
python -W error scalar_cold_packing_native_power_audit_reference.py "$NATIVE_OUT"
~~~

The comparisons are linear transfer controls at known
fixed inputs, not a Lyman-alpha or galaxy likelihood.

## Exact-mode comparison

| Dense arm at $z=5$ | $k=.3$ | $k=1$ | $k=2$ | $k=5$ |
|---|---:|---:|---:|---:|
| $F=1$: $P/P_{\rm CDM}$ | .98090 | .70259 | .25956 | .09922 |
| $F=200000$: $P/P_{\rm CDM}$ | .999994 | .999888 | .999408 | .994840 |

Wavenumbers here are in $h\,{\rm Mpc}^{-1}$.
Over the six retained modes and both redshifts,
the larger joint-tier relative drift is .763% for
$F=1$ and .0126% for $F=200000$. The matched
pure-CDM grids differ by less than $2.6\times10^{-9}$
relative at the same modes. The continuum number
normalizations match across the two shapes and tiers,
while the density-weighted relativistic moment drops:
the ideal early ultrarelativistic equivalent
$\Delta N_{\rm eff}$ is .01398 and .0002391 respectively
in the dense target shapes. These are input-shape
moments, not an inferred observational constraint.

[Exact modes and independent normalization checks](results/scalar-cold-packing-native-power-2026-11-26.json);
[compact raw power grids](results/scalar-cold-packing-native-compact-raw-2026-11-26.json.gz),
[dense raw power grids](results/scalar-cold-packing-native-dense-raw-2026-11-26.json.gz).
