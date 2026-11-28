# Independent bracket for finite normalized cold-profile separation

This check validates the exact-bin normalized-profile LP without using its nominal equality multipliers or objective value as a certificate. The LP supplies a candidate signed support $s$ with $|s_i|\le 1/2$ and two proposals for cell-simplex weights. For either chart, a lower bound on the normalized signed response comes from the independent support of its positive panel mixtures. For any trial $t$, let $h=\min(U-tT)$, where $U$ is the signed linear response and $T$ is the unsigned total. The minimum separates over the raw cell simplexes, and $T$ has a positive, independently computed lower bound $T_{min}$. Then $t+\min(0,h)/T_{min}$ is a lower bound on $\min U/T$. We apply this to $s$ on the first chart and $-s$ on the second. The sum lower-bounds the TV distance. Independently projected nonnegative weights and separately normalized profiles give an upper witness. Both include small floating-point guards. No directed rounding or continuous-momentum claim is made.

On their pair-specific grids, the six S2/S4 and S4/S8 observable checks pass both q meshes and have absolute brackets no wider than 7.09e−10 for either moment restriction. With each chart's own cold cell count and mean excess held, the minimum finite-profile TV is:

| Pair | Forward number | Inverse number | Forward excess |
| --- | ---: | ---: | ---: |
| S2/S4 | 2.0225559% | 2.0203947% | 1.8711101% |
| S4/S8 | 0.5377455% | 0.5329571% | 0.4942920% |

The accompanying 17-bin common S2/S4/S8 comparison has brackets no wider than 4.84e−10. These checks corroborate the existing finite-panel values via a different verification, rather than providing a new evolution, a spatial asymptotic error bound, or a physical distribution distance.

From `kinetics/`, for either adjacent accepted pair, with `a=2; b=4` or `a=4; b=8`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar; s=results/early-pole-fine-charge2-scalar
arc=results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
for ch in forward_number inverse_number forward_excess; do
 python -W error early_pole_endpoint_cold_profile_exact_finite_distance.py --independent-bracket \
   --levels "$a" "$b" --left "$p$a"-endpoint-continuum-cold-{1,2}.jsonl.gz \
   --right "$p$b"-endpoint-continuum-cold-{1,2}.jsonl.gz \
   --states "$s$a"-x05-state-"$d".json.gz "$s$b"-x05-state-"$d".json.gz --archive "$arc" --observable "$ch"
done
python -W error early_pole_endpoint_cold_common_grid_spatial_ladder.py --independent-bracket \
   --s2 "$p"2-endpoint-continuum-cold-{1,2}.jsonl.gz \
   --s4 "$p"4-endpoint-continuum-cold-{1,2}.jsonl.gz \
   --s8 "$p"8-endpoint-continuum-cold-{1,2}.jsonl.gz \
   --states "$s"{2,4,8}-x05-state-$d.json.gz --archive "$arc"
```
