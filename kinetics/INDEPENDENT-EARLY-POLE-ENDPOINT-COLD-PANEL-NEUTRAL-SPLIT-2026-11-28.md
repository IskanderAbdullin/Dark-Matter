# Does the conditional energy/number neutrality ordering survive finite reweighting?

Use the same S2 and S4 positive finite-panel deformations as the [count/mean excess study](INDEPENDENT-EARLY-POLE-ENDPOINT-COLD-PANEL-REWEIGHT-ENVELOPE-2026-11-28.md). Here optimize a bound on (n_E/n_N-1), *not* on the two neutral roots separately. Fix the charge distributions, scalar mass and photon pole.

Writing (B_N,D_N) for inverse and slope number moments and (B_X,D_X) for the corresponding excess-energy moments gives
`n_E/n_N - 1 = (B_X D_N - B_N D_X) / [B_N (m D_N + D_X)]`.
Normalize to the four original moments. Expand the numerator into its baseline, correlated linear change, and quadratic remainder. The linear part uses the same finite nonnegative count/mean polytope for all four components; then bound the quadratic remainder and positive denominator conservatively. Also rewrite the quadratic exactly into products with the paired small differences `delta X - delta Y` and `delta D - delta B` in normalized units. Optimize these two difference factors linearly on the same feasible cell vertices and intersect the two valid quadratic bounds. A synthetic all-vertex/interior fixture checks this algebra. Both original panel orders give consistent envelopes; neither numerical result is rounded in a directed way.

| Chart | Count-only allowed interval for (n_E/n_N-1) | Count and cell mean excess |
| --- | --- | --- |
| S2 | [+1.330e−8, +7.448e−8] | [+3.715e−8, +3.995e−8] |
| S4 | [+2.859e−8, +4.398e−8] | [+3.462e−8, +3.475e−8] |

The looser S2 remainder that boxes the two bulk products independently allows zero; it does not exhibit a sign-flipped configuration. Keeping the pair differences in the quadratic shows positive ordering on **both** fixed supports even in the count-only case. Adding the per-cell mean constraint narrows the two positive ranges substantially. This is stronger than separately bounding (n_E) and (n_N), which discards their correlation, but it cannot be used beyond the admitted panel values or under changes of optical parameters or actual evolved spectra. It does not predict a physical photon occupation.

Reproduction from `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_cold_panel_neutral_split_bound
d=2026-11-28
for s in 2 4; do
 python -W error early_pole_endpoint_cold_panel_neutral_split_bound.py \
  results/fine-charge2-scalar$s-endpoint-continuum-cold-{1,2}.jsonl.gz \
  results/early-pole-fine-charge2-scalar$s-x05-state-$d.json.gz \
  results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json --scalar-div "$s"
done
```
