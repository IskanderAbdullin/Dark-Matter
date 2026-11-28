# Raw CLASS momentum-node phase control

In the [pressure-balanced control](INDEPENDENT-CLASS-PRESSURE-BALANCED-CONTROL-2026-11-26.md),
smoothing and the number of nodes changed together. This test does
**not** smooth or alter the fine input distribution. Run the same
native synchronous CLASS setup with 128 truncated-trapezoid nodes
and change only the maximum sampled q from 20 to 21.3702247191.
Both cutoffs lie beyond the populated export's recommended q limit
(19.8367). The shifted endpoint was selected using a scan of
quadrature moments and dilation residuals, **not** the power ratio.

| Raw fine input, 128 nodes | q maximum 20 | q maximum 21.3702 |
| --- | ---: | ---: |
| Numerical physical number residual | +0.2093% | −0.07983% |
| Dilation identity residual for the number moment | +0.5412% | +0.01426% |
| Maximum absolute deviation from matched CDM, z=0 | 0.4521% | 0.01205% |
| Maximum absolute deviation from matched CDM, z=5 | 0.4521% | 0.01202% |

The only change in these two native inputs is qmax (plus the output
path). The maximum difference of their power ratios is 0.4401
percentage points across 95 comparison modes, k=0.01346–2.58357 h/Mpc.
At the unshifted endpoint the power excess also appears on a
denser k grid. In contrast, the earlier 64-node raw run shifted
the power downward by as much as 0.6100%. Thus these percent-level
features are numerical in these controls. A small residual on one
selected quadrature is **not** evidence for or a bound on the
physical transfer. It still lacks independent momentum convergence
and production-era perturbations.

Use the reproduction setup of the pressure-balanced test. In
addition to its original CDM and smooth reference cases, run
the generated configurations for psd_true_sync_fine,
psd_raw128_matched_sync and psd_raw128_q2137_matched_sync.
Then compare the requested files by name:

~~~sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for c in psd_true_sync_fine psd_raw128_matched_sync psd_raw128_q2137_matched_sync; do
  /path/to/class "native/$c.ini" >"native/$c.log" 2>&1
done
python -W error archived_class_pressure_balanced_reference.py native \
   --compare psd_true_sync_base psd_raw64_matched_sync psd_true_sync_fine \
             psd_raw128_matched_sync psd_raw128_q2137_matched_sync
~~~
