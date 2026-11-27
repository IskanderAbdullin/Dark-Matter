# Matching both species' retained number and bare energy leaves most of the finite contrast

Extend the [first-window hard-input swap](INDEPENDENT-EARLY-POLE-PREPARED-HARD-SWAP-TRAJECTORY-2026-11-27.md) at the same 16-step width and initial cold and photon data. Retune the continuous parent and daughter separately, using two common affine changes to the FD logits of each chosen species. Require each substituted spectrum to match its own native retained number and bare energy before evolution. Preserve the charged support and reinitialize the active response at that state. In the two-retuned arm the parent residuals are (0, −1.11e−16) and the daughter residuals are (+3.47e−18, 0), respectively.

| Prepared input changed from native | Starting cold rate | First-interval mean | Fraction of original finite gap retained |
| --- | ---: | ---: | ---: |
| Neither | .07863619025 | .06694243427 | 0 |
| Parent retuned, daughter native | .07283452112 | .06643447117 | 30.72% |
| Daughter retuned, parent native | .06302610960 | .06611548916 | 50.01% |
| Parent and daughter both retuned | .05722471602 | .06560747061 | 80.73% |
| Both continuous, no retuning | .05690395094 | .06528878623 | 100% |

Joint retuning leaves 98.52% of the instantaneous contrast and 80.73% of the finite contrast. Thus the retained number and bare energy of the two charged species do not by themselves determine these cold rates. The matched daughter result is not consistent with explaining its rate shift by its small original bare-energy difference alone. This is one manufactured retuning, not a unique split between shape and moments. It does not hold interaction energy, screening moments or the entire phase-space density fixed. All matched trajectories pass the unchanged terminal conservation checks.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
for k in matched-daughter matched-both; do
  python -W error early_pole_prepared_hard_swap_trajectory_reference.py "$a" --case "$k" --steps 16
done
```

Raw: `results/early-pole-prepared-hard-joint-moment-match-2026-11-27.json.gz`.
