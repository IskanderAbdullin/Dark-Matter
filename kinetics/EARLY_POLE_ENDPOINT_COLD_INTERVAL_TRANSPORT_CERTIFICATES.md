# Positive but small interval-transport distance across the cold spatial ladder

A normalized cold profile on one chart assigns mass to its native q panels; it does not locate the response within a panel. Match these positive masses to the other chart, charging only the distance between the two intervals. This can cost zero for overlapping or touching supports. Compare unchanged panel totals, reweighting with count fixed, and reweighting with both count and mean excess fixed independently in each chart's own cells.

A linear optimizer proposes a coupling and potentials. For an independent upper bound we project to positive shell simplexes, renormalize and repair the coupling to its two resulting marginals. For a lower bound we correct the potentials so their sum never exceeds a local transport cost, then separately minimize their weighted responses over each chart's normalized finite family. This gives a checked bracket with ordinary floating-point guards, not directed rounding. The original unscaled nominal path failed its strict residual gate on several cases; no results are taken from those failed proposals.

On both q meshes and all three channels, count-only lower bounds are zero, with checked upper witnesses below 2.33e−15 in q for this interval problem. With each chart's own count and mean excess fixed the costs are:

| Adjacent levels | Forward number | Inverse number | Forward excess |
| --- | ---: | ---: | ---: |
| S2/S4 | [7.9375,7.9377]e−11 | [8.4455,8.4457]e−11 | [2.4327,2.4328]e−10 |
| S4/S8 | [6.4327,6.4352]e−12 | [6.5926,6.5951]e−12 | [1.1467,1.1470]e−11 |
| S8/S16 | [0.9685,1.2100]e−14 | [0.9903,1.2393]e−14 | [2.2066,2.4392]e−14 |

The last pair retains a resolved positive cost after count and mean fixing, despite the arbitrary within-panel positions. It is in q units, whereas the fixed-bin TV statistic is dimensionless and penalizes redistributions that may require almost no displacement. Neither statistic establishes a physical spectrum, a true momentum-transfer rate, or a spatial convergence exponent.

From `kinetics/`, set a to 2, 4 or 8 and b to twice a:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
for c in forward_number inverse_number forward_excess; do
 python -W error early_pole_endpoint_cold_profile_interval_transport.py --independent-bracket \
  --levels "$a" "$b" --observable "$c" \
  --left "$p$a"-endpoint-continuum-cold-{1,2}.jsonl.gz \
  --right "$p$b"-endpoint-continuum-cold-{1,2}.jsonl.gz \
  --states "$s$a"-x05-state-$d.json.gz "$s$b"-x05-state-$d.json.gz \
  --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
done
```
