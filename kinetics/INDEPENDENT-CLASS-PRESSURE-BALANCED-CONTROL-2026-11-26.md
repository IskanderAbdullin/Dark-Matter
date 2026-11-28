# CLASS: pressure-balanced control for the archived cold stable component

This is a numerical sensitivity test in synchronous gauge for the
**retained 20 MeV stable fermion** from the low-reheat paired archive.
It omits the scalar population and exterior cohorts. Starting
perturbations are standard adiabatic perturbations of an extrapolated
collisionless spectrum. The calculation does not evolve perturbations
during production. A tiny `ω_cdm=10⁻⁸` anchors the synchronous gauge;
the remaining `ω=0.11816900374571636` is fixed across all finite
momentum tests. Matched CDM uses the same total.

## Change the input without changing number or the late pressure moment

Smooth the log-number measure $q^3 f_0(q)$ by a Gaussian in
$\ln q$ with width 0.055. Shift $\ln q$ down by $0.055^2$
and normalize back to the same number. In the continuum this keeps
$\langle q^2\rangle$ invariant: the Gaussian multiplies the
unshifted $\langle q^j\rangle$ by $\exp(j^2\sigma^2/2)$.
Here the normalized moments of order 1 through 4 change by
−0.151%, roundoff, +0.455%, and +1.217%. This preserves the
leading nonrelativistic pressure moment; it does not preserve the
distribution or exact relativistic pressure. In a direct
background integral from (a=10^{-13}) to 1, the largest
fractional changes in the species' density and pressure are
0.151% and 0.187%. The density change is about 4.1e-8
at the first sampled scale factor after the archived snapshot.

At 256 manual momentum nodes the unbalanced fine export
has relative number quadrature error 7.42e-4 and
maximum dilation-identity residual 0.0125. The modified export
gives 7.02e-5 and 4.00e-5, respectively.
These are direct quadrature diagnostics against each spline's own
continuum moments, not error bars on CLASS power.

## Native result

Use unmodified CLASS v3.3.4. Compare 95 near-matched (k) samples
between 0.01346 and 2.58357 h/Mpc.
Maximum $|P/P_{\rm CDM}-1|$ over those modes:

| Tabulated input and settings | z=0 | z=5 |
| --- | ---: | ---: |
| Standard narrow export, 64 q nodes, original k sampling | 0.0032792 | 0.0032794 |
| Fine narrow export, 64 q nodes, matched k sampling and fluid trigger | 0.0060999 | 0.0060994 |
| Pressure-balanced fine export, 256 q nodes, matched k sampling and fluid trigger | 0.00003490 | 0.00003464 |

The two matched CDM/pressure-balanced (k) grids differ by at
most 1.49e-8 fractionally. Interpolating the balanced
case across that tiny shift changes its power by at most
3.51e-8 fractionally. Replacing the matched CDM
reference with an interpolated, differently sampled CDM run instead
moves the reference by as much as 1.85e-4.

The raw finite-grid deviations are not robust across these controls.
Do not interpret the smoother's smaller residual as a bound on
the undistorted physical transfer. These cases vary momentum
resolution and smoothing together. The exactly matched raw 256-node
comparison, fluid-hierarchy convergence, and production-era
perturbations remain separate checks.

## Reproduce

Build unmodified [CLASS v3.3.4](https://github.com/lesgourg/class_public/tree/v3.3.4)
and use the [official radiation table as described in the setup](INDEPENDENT-BENCHMARK-SETUP.md).
Work from this directory:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
snap=results/low-reheat-quantum-bridge-wide-shape-2026-11-26.json
python -W error export_archived_class_psd_reference.py EOS.dat "$snap" standard --sigma .010 --nodes 3200
python -W error export_archived_class_psd_reference.py EOS.dat "$snap" fine --sigma .006 --nodes 6200
python -W error archived_class_pressure_balanced_reference.py native \
    --standard-export standard --fine-export fine --audit input-audit.json
for c in cdm_true_sync cdm_balanced_grid_sync psd_true_sync_base psd_raw64_matched_sync psd_balanced_sync; do
  /path/to/class "native/$c.ini" >"native/$c.log" 2>&1
done
python -W error archived_class_pressure_balanced_reference.py native \
    --compare psd_true_sync_base psd_raw64_matched_sync >power-audit.json
```
