# Fixed bare daughter moments leave the initial suppression but weaken the finite rebound

Keep the same archived held-response graph and the two substituted
local daughter cells exactly unchanged. Use the
[previous two-moment construction](INDEPENDENT-EARLY-POLE-DAUGHTER-COMPLEMENT-RETUNE-2026-11-27.md)
to compensate the retained daughter number and bare energy in the
other 18 cells. The evolved construction reproduces the static
fit and all recorded preparation residuals. The cold seed, parent
and photon inputs and hard support remain fixed.

| Local prepared input | Initial shift from native | First-window mean shift from native |
|---|---:|---:|
| Unmatched | -0.01551056101 | +5.17245309e-4 |
| Bare number and energy matched | -0.01551107287 | +8.85152056e-5 |

The matched variant alters the initial shift by only about 33 ppm
but reduces the finite positive effect by 82.9%. It stays positive
in this run. Energy, entropy and total hard number pass the original
unsplit gates. The cumulative difference from the unmatched path is
negative at every saved positive time, reaching -1.4291e-10 in cold
count. This establishes sensitivity of the finite effect to this
specific outside compensation, not a unique causal split between
bare moments and shape. The compensating operation changes other
moments and occupations as well. A complementary outside-only
control under the same logit change can test finite additivity.

Reproduce with `early_pole_joined_simultaneous_chart_reference.py`
and the paired archive plus `--retune-local-daughter`, then pass the
native, unmatched local, matched local and static retune receipts
to `early_pole_local_daughter_retuned_time_compare.py`.
