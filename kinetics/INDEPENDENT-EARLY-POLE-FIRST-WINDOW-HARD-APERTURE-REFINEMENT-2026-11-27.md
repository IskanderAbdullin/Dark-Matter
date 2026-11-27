# Refining hard cells changes the early source even when retained moments agree

Insert four common anchors in each charged hard grid near the original aperture. The charged endpoints, photon grid, cold scalar seed and overall retained species counts are unchanged. Compare 20 to 24 cells per charged species for the original first interval with 16 time steps and eight cold bins, preserving the original five channels and solver gates. Separately evaluate the initial source with aperture subdivision factors 1, 2, 4 and 8; these static deeper calculations are not later-time simulations.

| Hard aperture subdivision | Native initial cold rate | Continuous initial cold rate |
| ---: | ---: | ---: |
| 1 | .07863619025 | .05690395094 |
| 2 | .07863564548 | .06054576741 |
| 4 | .07863542794 | .06032821588 |
| 8 | .07863542212 | .06034537644 |

The static checks keep retained species totals and bare energies within 1.2e−16 across the subdivisions. The original continuous initial rate is therefore sensitive to how those same moments are resolved within the aperture, increasing about 6.4% at factor 2 and then changing less at the deeper factors. The native starting rate is nearly unchanged. At finite exposure the tested 1→2 refinement increases the mean rate by 575 ppm in native and 2,017 ppm in continuous. The first-mean contrast relative to native is 2.470% at 20 cells and 2.330% at 24 cells. The accepted terminal defects still satisfy the fixed energy and entropy thresholds. Neither the static deeper sequence nor this one finite spatial refinement establishes a continuum rate.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
for p in steps cell-conservative; do
  for n in 1 2; do
    python -W error early_pole_cold_screened_split_reference.py fine \
      --duration 3.3333333333333335e-7 --count 16 --pair --same-soft-strength 1 --compton-strength 1 \
      --conversion-strength 1 --screening-update full --projection "$p" --snapshots \
      --elastic-energy-budget 5e-15 --aperture-subdivision "$n" > "$p-hard-sub$n.json"
  done
  python -W error early_pole_first_window_hard_aperture_refinement_reference.py "$a" \
    "$p-hard-sub1.json" "$p-hard-sub2.json" --projection "$p"
done
```

Raw: `results/early-pole-first-window-hard-aperture-refinement-2026-11-27.json.gz`.
