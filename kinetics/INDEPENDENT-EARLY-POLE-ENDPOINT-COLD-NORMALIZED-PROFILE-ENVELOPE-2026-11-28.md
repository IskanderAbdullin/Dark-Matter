# A cold-profile difference survives independent positive reweights on both charts

The matched S2/S4 gross forward profiles differ in normalized mass on their [exact shared q partition](INDEPENDENT-EARLY-POLE-ENDPOINT-COLD-COMMON-Q-PARTITION-2026-11-28.md). Test whether those differences could disappear by independently reshaping each original scalar shell at its own held pole. Use the finite positive per-panel occupations and the known cell amounts from the [finite reweight experiment](INDEPENDENT-EARLY-POLE-ENDPOINT-COLD-PANEL-REWEIGHT-ENVELOPE-2026-11-28.md). Allow any nonnegative scalar panel weights preserving each in-band cell count, then also each in-band cell mean excess. No interpolation, charge-shape change or altered optical input within an arm.

For each exact common group, bound its normalized gross-forward share by optimizing group minus (r) times the **same** total on each positive vertex family and solving for the extremal ratio. Do not divide independent extrema of the group and total. If the share intervals for the two independent charts are disjoint in a group, the gap gives a verified sign and a minimum gain or loss. The maximum of the two summed one-sided gaps is a valid lower bound on normalized total variation, as all gains must balance losses. It need not be the exact optimal jointly attained distance.

- The original shared-bin profile difference is 2.02402% total variation. Allowing independent count-preserving shapes leaves a **0.72914%** lower bound. With both per-cell count and mean excess fixed, the lower bound is **2.02255%**.
- In the key bin [3.30599941e−6,3.97454049e−6], the S4 loss versus S2 is at least 0.64936 percentage points even with counts alone; with means fixed at both references it is at least 2.02255 points.
- Exact support, panel weighting and the two independent collision quadratures are checked. These are finite, frozen, conditional comparisons. Neither yield agreement nor these moment constraints establish a physical cold spectrum or a true continuum numerical error bar.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_cold_normalized_profile_envelope
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
python -W error early_pole_endpoint_cold_normalized_profile_envelope.py \
 --left "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --right "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$s"{2,4}-x05-state-$d.json.gz --levels 2 4 \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
```
