# The initial directional signs survive packet refinement, but one magnitude does not

Repeat the [prepared-state initial directional audit](INDEPENDENT-EARLY-POLE-INITIAL-COLD-RATE-DIRECTION-2026-11-27.md) with different elastic packet quadrature orders while holding the hard cells, the scalar seed, and all other channels fixed. The tuple is (radial, transfer, angular) in the graph constructors. Report the center probe (3e−15) at each setting; the nearby signed probes and additivity checks are retained (center-mode relative additivity defect below 7e−7). These are local derivatives of the cold source, not rates averaged over a finite interval.

| Unlike and same quadrature | Native combined derivative | Continuous combined derivative |
| --- | ---: | ---: |
| (1,1,2), original | −3.5570e6 | +2.0431e4 |
| (1,1,4), angular refinement | −3.4200e6 | +2.0697e4 |
| (2,1,2), radial refinement | −3.8673e6 | +9.9700e4 |
| (2,1,4), both changes | −3.7784e6 | +9.1707e4 |

The signs agree across these tests, but the continuous magnitude changes by almost a factor of five. The unlike and same directions partially cancel there; refining only the unlike radial quadrature moves its total derivative to +9.4280e4, versus +2.5850e4 if only the same-flavor quadrature is refined. The unlike change accounts for most of the observed variation in the signed sum. The initial derivative magnitude is not established to be converged, and this test alone does not determine the corresponding finite-time gain change. Default probe output has been replayed byte for byte.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
for p in steps cell-conservative; do
  for orders in '1 1 2' '1 1 4' '2 1 2' '2 1 4'; do
    python -W error early_pole_joined_initial_rate_direction_reference.py "$a" --projection "$p" \
      --unlike-quad $orders --same-quad $orders
  done
done
python -W error early_pole_joined_initial_rate_direction_reference.py "$a" \
  --projection cell-conservative --unlike-quad 2 1 2
python -W error early_pole_joined_initial_rate_direction_reference.py "$a" \
  --projection cell-conservative --same-quad 2 1 2
```

Raw: `results/early-pole-first-slope-elastic-packet-refinement-2026-11-27.json`.
