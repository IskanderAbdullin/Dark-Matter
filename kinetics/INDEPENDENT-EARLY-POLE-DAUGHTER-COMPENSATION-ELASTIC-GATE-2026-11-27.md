# Removing elastic scattering restores the initial local suppression

Compare the same native, local, complement-only and combined
preparations on the simultaneous chart, once with its two elastic
terms and once with those terms omitted. Hold the common exposure,
selected hard supports, seed, tolerance and all original
preparation data identical. In both runs retain reciprocal decay,
hard pair, Compton and conversion. Modeled dressing remains; this
is not a change of the imposed dielectric table or a comparison
of total rates across different physical models.

All reported shifts are relative to the native calculation
with the same channels. The solver's original gates apply to all
eight paths. The complement metadata is checked byte-equivalent
as a parsed payload across the channel switch.

| Prepared change | Joined mean shift | Nonelastic mean shift |
|---|---:|---:|
| Local only | +5.17245e-4 | -1.55091e-2 |
| Complement compensation only | -4.28893e-4 | -4.20643e-7 |
| Local plus compensation | +8.85152e-5 | -1.55096e-2 |

Without elastic scattering the local effect remains close to the
initial suppression (-1.55106e-2) instead of turning positive.
The complement-only shift is just 0.0981% as large as its joined
value, and adding the complement on top of local gives 0.1227%
of the corresponding joined difference. The absence of two
elastic terms thus changes the finite response by much more than
the small static tail differences. The comparison is at the same
dielectric graph and original short exposure.

The omitted-channel comparison is a conditional numerical test.
It does not identify a unique contribution per elastic reaction or
estimate the cosmic production of a physical relic.

From `kinetics/`, repeat the previous four simultaneous controls
adding `--without-elastic`, then feed the same four nominal
prepared data into
`early_pole_daughter_compensation_elastic_gate_compare.py`
under each option, in native, local, complement-only and combined
order.
