# The held composition and the absolute direct abundance are different normalizations

The [previous few-ppm matter-power checks](INDEPENDENT-WARM-DIRECT-NATIVE-WIDE-CONVERGENCE-2026-11-26.md)
hold the ratio of the independently produced scalar count
to the archived stable count, then set the **total** dark
density to omega_dark=.12. That is a valid comparison of
a deliberately rescaled dark composition. It is **not**
the fraction obtained by keeping the absolute direct
scalar yield and independently providing the remaining
omega_dark with CDM. Separate those two scenarios
before interpreting the native effect.

| Conditional normalization | Result |
| --- | ---: |
| Held q64 stable phase, if it simply redshifts without count changes | omega_stable=75.7317 |
| Held q64 scalar phase under the same continuation | omega_scalar=.000199355 |
| Numerical direct lepton plus measured pion vacuum yield | Y_scalar=6.69938e−6 |
| Absolute omega_scalar of that independent direct yield | .000158975 |
| Ratio of direct scalar to direct scalar plus held stable | 2.099185e−6 |
| Fraction if the absolute direct yield instead replaces part of omega_dark=.12 | .001324794 |

The first row is a deliberately conditional extrapolation
of an overabundant archive, not a prediction or a viable
stable cosmology. Renormalizing the whole composition in
the earlier scenario reduces absolute scalar and stable
amounts by a common factor of about 631.1. In the second
scenario the stable CDM share is prescribed independently
to fill the remainder; it is **not** predicted by the same
portal model. Thus neither case proves a fully solved
production history. The conversion between Y and omega
uses the same present entropy density; a separate SI
calculation verifies the unit conversion.

## Absolute envelope for the restricted direct channels

Use a positive invariant-mass envelope to bound the
full direct scalar yield from 50 MeV to .1 MeV, no
cold cut. Integrate the MB total lab energy exactly
through the Bessel K1 identity. FD is bounded above by
MB for charged leptons, and the physical pion mass gives
a finite uniform Bose multiplier at each temperature.
Use upper source and lower H/entropy endpoints on
independently refined grids. Include the measured KLOE
fixed-bin choice plus twice the quoted marginal errors,
and the two independent low positive moment bounds.
For the still missing charged-pion high band use
the explicitly assumed R_pi≤1, not an extrapolated
measurement. The direct scalar upper omega is
.00022551 under these conditions, or .0018793 of
omega_dark=.12. Other sources and dark conversion
or rescattering are not included.

For scale, assigning half of the reference dark density
to the same spectrum would be about 377 times its
numerical absolute direct yield. This is not an
inferred coupling rescaling or production pathway.
Native inputs for the absolute-yield scenario and
the deliberately rescaled scenarios are kept separate.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_direct_vacuum_absolute_source_envelope_reference.py EOS.dat KLOExSec.data \
 "$D/warm-late-cold-pion-dispersive-recoil-2026-11-26.json" "$D/warm-direct-vacuum-scalar-psd-2026-11-26.json" > absolute.json
python -W error warm_direct_vacuum_native_absolute_abundance_reference.py EOS.dat \
 "$D/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz" absolute.json higherk absolute
~~~

[Absolute vacuum production envelope](results/warm-direct-vacuum-absolute-source-envelope-2026-11-26.json);
[reference abundance scenarios](results/warm-direct-vacuum-two-normalizations-2026-11-26.json).
