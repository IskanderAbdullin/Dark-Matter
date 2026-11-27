# Check which parts of the fixed-metric response actually change

This follows the [positive-shell fixed metric fold](INDEPENDENT-WARM-DIRECT-VACUUM-FIXED-METRIC-2026-11-26.md).
No new gravitational metric is generated.

## Change the production clock without changing its peak

Build a full positive shell spectrum for each of the fixed
50 MeV clock controls. The radiation reference gives
6.69938e−6 total scalar yield in the two selected vacuum
sources. The alternate duration factors at 30 and 20 MeV
give .58857 and .43924 of this yield. Keeping only 20%
of the original proper time above 45 MeV gives .68933 of
the total yield, but only .31548 of the p/T<.1 yield at
equality. The largest change across the normalized CDF
in that last control is .01933.

The same production kernels and entropy-to-momentum map
apply throughout. Keep the recorded metric fixed to
separate the spectrum effect from overall scalar number.
At k=6.736 Mpc⁻¹ the last clock control changes the
normalized scalar shape amplitude from .46303 to .45327.
At fixed scalar count this is −2.05e−8 in a CDM-plus-scalar
matter-amplitude proxy. Letting the scalar count change
as well changes the sign and gives +3.36e−7, because
less of the mixture is assigned the more suppressed
scalar. This is a proxy with the stable component held,
not a calculation of a different cosmology.

## Check driver choice and the low-q part of the response

Repeating the same scalar shapes on the alternative
stored metric fixtures with the earlier fluid grouping
enabled changes the individual ratios by less than
2.01e−8 in the two fixtures. This controls that narrow
driver replacement, not general metric uncertainty.

Separately remove the q/T0≤.1 band from each shape and
renormalize its complement. The removed energy fraction
at the final time is about 1.867e−4 for the direct mix
and 4.976e−4 for the held shape. The low-q band alone
responds more like CDM on these two modes, but the
difference between the whole direct and held scalar
amplitudes changes by at most 1.98e−4 after removing it.
The .05835 and .07627 shape contrasts are mainly carried
by the broader distributions, not the rare equality-soft
tail. This is an algebraic decomposition of the same
linear functional, not a solved initial condition for
independently clipped physical particles.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_direct_vacuum_clock_shape_response_reference.py EOS.dat KLOExSec.data \
 "$D/warm-direct-vacuum-scalar-psd-2026-11-26.json" \
 "$D/warm-late-cold-source-temperature-localization-2026-11-26.json" \
 "$D/warm-direct-vacuum-fixed-metric-response-2026-11-26.json" "$D"/los-metric-off/*.npz
python -W error warm_direct_vacuum_response_band_audit_reference.py EOS.dat \
 "$D/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz" \
 "$D/warm-direct-vacuum-scalar-psd-2026-11-26.json" "$D"/los-metric-off/*.npz
# Repeat the fixed-metric response doc's command using "$D"/los-metric-on/*.npz.
~~~

[Clock/shape decomposition](results/warm-direct-vacuum-clock-shape-response-2026-11-26.json);
[low-q decomposition](results/warm-direct-vacuum-response-band-audit-2026-11-26.json);
[alternate metric driver](results/warm-direct-vacuum-fixed-metric-on-2026-11-26.json).
