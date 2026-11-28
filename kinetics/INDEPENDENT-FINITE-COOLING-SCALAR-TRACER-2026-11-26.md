# Native scalar tracer for the finite cooling-history prompt source

This tests the two [resolved distributions](INDEPENDENT-FINITE-COOLING-PHASE-SPACE-2026-11-26.md)
at actual prompt-only scalar densities, in unmodified CLASS 3.3.4.
The 86 eV scalar is evolved collisionlessly and adiabatically; the
20 MeV dominant species is represented as ordinary cold dark matter
with Ωh²=0.120. The scalar densities at the nominal rest count are
2.36354×10⁻¹² (instantaneous) and 2.23256×10⁻⁹ (finite with Tmax/TR=5).
This is a comparison of two **scalar tracer density contrasts**,
not a comparison of physically dominant scalar transfer functions.

Two positive reconstructions, log widths 0.10 and 0.14, preserve
the recorded number and squared momentum. They differ in higher
moments and use 160 and 224 native indefinite momentum nodes.
The colder tracer requires coordinate scales 0.40277566 and
0.33892848, chosen by the independent moment/dilation screen;
the physical distribution is invariant under the coordinate
choice. The corresponding screened maximum defects are 3.43×10⁻⁷
and 4.41×10⁻¹⁰, while the unremapped native instantaneous variants
pass their moment gate on the original scale. The solver's own
density agrees with the scalar input; the largest k times native
conformal time at the snapshot is 5.57×10⁻⁴, and all included modes
begin their native perturbations later. This is a check on the
native clock, not matching through the actual reheating era.

| Redshift | Instantaneous scalar/CDM at k=10 h/Mpc | Finite cooling scalar/CDM at k=10 h/Mpc | Instantaneous first amplitude half-crossing |
| --- | ---: | ---: | ---: |
| 0 | 0.60658 | 0.97865 | not reached |
| 5 | 0.33479 | 0.93278 | 6.1818 h/Mpc |
| 20 | 0.15897 | 0.84709 | 3.2640 h/Mpc |

No finite-cooling amplitude half-crossing is found in the simulated
range (comparison reported up to 10 h/Mpc). Across 0.015–10 h/Mpc
and the three redshifts, the two shape/node controls differ
absolutely by at most 1.07×10⁻⁴ for instantaneous and 1.90×10⁻⁵
for finite. Both use the same perturbation hierarchy
(lmax 14 and ncdm fluid trigger 20); these are not hierarchy
convergence bounds. The separate late forward-pair diagnostic
shows very small absolute counts but a possible comparable
addition to the instantaneous **prompt scalar** count, which is
not included in these spectra. No power suppression of total
matter, early entropy or isocurvature perturbations, or observational
constraint follows from the tracer comparison.

## Native hierarchy extension

With identical input, clock and momentum nodes at σ=0.10, also run
lmax=22 and fluid trigger 42. The maximum absolute difference
in the scalar/CDM contrast on the common k interval is
1.70×10⁻³ (instantaneous) or 1.20×10⁻³ (finite), with the
largest difference at z=20. At k=10 and z=20, the finite case
changes by −1.06×10⁻⁴; the instantaneous case changes by
4.1×10⁻⁷. These are comparisons of two closure choices, not a
no-fluid convergence envelope. The primary table remains the
original hierarchy for reproducibility.

## Reproduction

Regenerate the PSD tables with the documented exporter, then use:

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error prepare_finite_cooling_scalar_tracer_reference.py     native --psd psd --manifest results/finite-cooling-class-psd-export-2026-11-26.json

for ini in native/finite_cooling_tracer_*_sync.ini; do
  OMP_NUM_THREADS=1 /path/to/class "$ini" >"${ini%.ini}.log" 2>&1
done

OPENBLAS_NUM_THREADS=1 python -W error prepare_finite_cooling_scalar_tracer_reference.py     native --eos EOS.dat --compare finite_cooling_tracer_instant_s10_sync     finite_cooling_tracer_instant_s14_sync finite_cooling_tracer_ratio5_s10_sync     finite_cooling_tracer_ratio5_s14_sync
~~~

For the hierarchy-only check, copy the two s10 INI files under
the matching names ending in s10_hierarchy_sync.ini. Change only the
root output prefix, l_max_ncdm to 22 and ncdm_fluid_trigger_tau_over_tau_k
to 42. Run both and compare with:

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_scalar_hierarchy_audit_reference.py native
~~~
