# Linear 86 eV scalar tracer in the archived low-reheat branch

Use the retained scalar distribution from the low-temperature
paired history. Its fixed-count density is
$\omega_\phi=2.98982\times10^{-10}$, about 2.53e-9 of the dark
density in this run. Treat the dominant stable population as CDM,
which avoids the [known stable-species numerical normalization issue](INDEPENDENT-CLASS-DILATION-INTERVENTION-2026-11-26.md)
in this independent scalar test. The distinction matters: this
does **not** compare a universe made of 86 eV dark matter to CDM.

Use unmodified CLASS v3.3.4, synchronous gauge, adiabatic initial
conditions for the late shape extrapolated backward without
production. Smooth the number measure in log-momentum by width
0.055 and shift log-momentum to preserve the leading $q^2$
pressure moment as in the earlier numerical control. The narrow and
first wide runs share the same scalar shape, 128 momentum nodes
and hierarchy. An additional wide run increases the multipole
cutoff from 10 to 14 and delays the fluid approximation
(criterion 12 to 20). Densities are fixed,
not normalized to a value that makes the scalar important.

| z | $\delta_\phi/\delta_{\rm CDM}$ at k=1 h/Mpc | at k=4 h/Mpc | at k=10 h/Mpc | First bracketed amplitude half, h/Mpc |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 0.970 | 0.817 | 0.562 | 11.99 |
| 5 | 0.909 | 0.602 | 0.292 | 5.41 |
| 20 | 0.802 | 0.384 | 0.132 | 2.86 |

The ratio uses both species' density contrasts from the same run.
The narrower and wider k calculations differ by at most 1.1e-4
in absolute ratio at their six overlapping probes. This is a
sampling control. With the later fluid switch, the maximum change
across their native wide grid is 7.87e-4 in absolute ratio, and the
three wide half-amplitude locations change by less than 0.0002 h/Mpc.
An independent raw-PSD control uses 160 nodes, qmax 63.8377
selected by external moments, and the later fluid switch.
It changes the wide ratios at k = 2, 4, 8 and 10 h/Mpc by
at most 9.9e-5 in absolute value relative to the balanced
late-switch run. Its tracer amplitude-half scales are
11.9890, 5.40664 and 2.85501 h/Mpc. This is evidence of
numerical robustness for the shown probes, not a general
proof of momentum or early-history convergence.
A higher-node narrow test compares raw 384, shifted raw 160,
balanced 256 and balanced 128 with a common k window. Across
the six probes the largest absolute ratio spreads are
8.0e-5, 3.8e-4 and 1.31e-3 at z=0, 5 and 20, respectively.
At k=1 h/Mpc the raw 384 values are 0.969771, 0.908875 and
0.801930. These extra controls do not extend the 384-node
result beyond k=1 h/Mpc.
At the archived scale factor,
the largest wide mode has kτ ≈ 0.00206 in this background; standard
adiabatic perturbations are initialized later with kτ ≤ 0.0701.
The scalar tracer increasingly follows CDM at late
times in this conditional linear history, but still loses
contrast on small scales. The amplitude-half crossing applies
only to this tracer ratio, not to matter power, an observed
structure scale, or a dominant scalar component. In linear
density weighting, the retained scalar is far too rare to
interpret this attenuation as a bulk matter-power suppression.

## Reproduce

Export the fine PSD from the existing low-reheat paired checkpoint
as in the pressure-balanced CLASS control, then:

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error archived_class_scalar_tracer_reference.py native --fine fine >input-audit.json
for c in scalar_tracer_smooth128_sync scalar_tracer_smooth128_wide_sync scalar_tracer_smooth128_wide_delayed_sync scalar_tracer_raw160_q638_wide_sync; do
  /path/to/class "native/$c.ini" >"native/$c.log" 2>&1
done
python -W error archived_class_scalar_tracer_reference.py native --fine fine \
  --compare scalar_tracer_smooth128_sync scalar_tracer_smooth128_wide_sync scalar_tracer_smooth128_wide_delayed_sync scalar_tracer_raw160_q638_wide_sync
~~~

To reproduce the higher-node narrow probe audit, run the generated
scalar_tracer_raw384_sync and scalar_tracer_smooth256_sync configs,
plus the existing 128 and shifted 160 narrow configs. Then use
\`\`\`sh
python -W error archived_class_scalar_tracer_reference.py native --fine fine --probe-only \
  --compare scalar_tracer_smooth128_sync scalar_tracer_smooth256_sync \
            scalar_tracer_raw160_q638_sync scalar_tracer_raw384_sync
\`\`\`
