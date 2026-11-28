# Transfer the unmatched charged-pion windows into the cold yield

The [KLOE recoil fold](INDEPENDENT-WARM-LATE-COLD-VISIBLE-KLOE-RECOIL-2026-11-26.md)
uses measured bins only on .10≤s≤.95 GeV². For the two
distinct regions below those bins, apply positive moments
already documented for the [theory threshold slice](INDEPENDENT-WARM-DISPERSIVE-THEORY-PION-THRESHOLD-2026-11-26.md)
and the [BABAR gap](INDEPENDENT-WARM-BABAR-DISPERSIVE-GAP-2026-11-26.md).
Use these only as assumed moment allowances on the same
vacuum bare convention, not as new point spectra or
confidence limits.

For each s box enclose the time-integrated kinematic
cold weight W(s) using a parent-angle box before the
Bose population factor. Cold eligibility requires
energy above the exact recoil gate. Pion BE product
is bounded from the minimum individual pion energy
and sum energy, before spin weights. The far total-energy
tail uses the conservative exp(−Q/T) integral times Q.
Divide the transfer maximum by a positive lower for
the muon dispersive weight H(s), and multiply by the
assumed whole positive moment; this needs no shape for
sigma(s) inside the two gaps. Use the full monotone EOS
for the opposite time endpoints from 50 to .1 MeV.

| Distinct interval in sqrt(s) | Assumed positive moment | Refined upper / held cold yield |
| --- | ---: | ---: |
| 2m_pi–.30 GeV (theory prior) | 1.3443e−8 | .005058 |
| .30–sqrt(.10) GeV (BABAR prior) | 6.69e−9 | .003816 |

The two independent allocations sum to less than 0.89%
of the held cold yield. They deliberately grant a
whole moment to each target slice and are conservative
even when not using information about overlap.

## Alternative: budget the overlapping measurements before allocating

If both whole moment caps and a chosen piecewise
KLOE **lower** spectrum apply simultaneously to a common
vacuum spectrum, first subtract the known positive KLOE
overlap. Use KLOE bins clipped at .50 GeV for BABAR, and
at .63 GeV for the theory prior. Pool what remains of the
theory moment between the two missing bands, and allocate
as much of it as BABAR allows to the band with larger W/H.
With a nonnegative lower spectrum equal to KLOE central
minus two quoted marginal errors per bin, the joint upper
becomes 4.788e−12, or 0.115% of the held cold yield.
Central-minus-one and central-minus-three checks give
0.0649% and 0.164%, respectively. This is a deterministic
overlap sensitivity: the binwise lower choices are not
proved statistical coverage or a fit to correlated data.

## Upper bands remain a separate conditional choice

For s≥.95 GeV², do not extrapolate the KLOE spectral
shape. A uniform assumed cap R_pi(s)≤Rmax gives a
conditional yield upper of 7.487e−4 Rmax of the held yield,
including the separate total-energy tail beyond sqrt(s)=4
GeV. With the lepton enclosure, the KLOE fixed-bin
central-plus-two-marginal choice and the two separately
allocated low bands, this gives a benchmark upper below
0.379 of the held yield at Rmax=1. A uniform cap of
Rmax≤830 would suffice for the conditional total of these
channels to stay below the held yield. Neither value is
inferred from the measured data, and there is no claim
about unmeasured other hadronic channels, in-medium effects,
a polarized parent or later rescattering.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
prefix="$D/warm-late-cold"
python -W error warm_late_cold_pion_dispersive_recoil_bound_reference.py EOS.dat \
 "$D/warm-babar-dispersive-gap-2026-11-26.json" \
 "$D/warm-dispersive-theory-pion-threshold-2026-11-26.json" \
 "$prefix-visible-exact-recoil-2026-11-26.json" > gap.json
python -W error warm_late_cold_pion_moment_overlap_reference.py gap.json KLOExSec.data > overlap.json
python -W error warm_late_cold_pion_high_band_headroom_reference.py EOS.dat \
 "$prefix-visible-exact-recoil-2026-11-26.json" "$prefix-visible-recoil-box-2026-11-26.json" \
 "$prefix-visible-kloe-recoil-2026-11-26.json" gap.json overlap.json
~~~

[Missing low bands](results/warm-late-cold-pion-dispersive-recoil-2026-11-26.json);
[optional overlap allocation](results/warm-late-cold-pion-moment-overlap-2026-11-26.json);
[high-band sensitivity](results/warm-late-cold-pion-high-band-2026-11-26.json).
