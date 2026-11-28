# Initially active cells do not carry the full finite-time daughter effect

Exchange continuous hard moments only in the [two cells per species](INDEPENDENT-EARLY-POLE-INITIAL-COLD-SOURCE-HARD-CELLS-2026-11-27.md) that contribute to the initial cold occupation contrast, leaving the rest of the native hard spectrum unchanged. Evolve the same five-channel first interval at 16 steps with the same conservative acceptance. The support, initial scalar and photon states, and active response initialization rules match the [full-region swaps](INDEPENDENT-EARLY-POLE-PREPARED-HARD-SWAP-TRAJECTORY-2026-11-27.md). Changing hard moments here recomputes the dressed geometry, unlike a static held-dispersion substitution.

| Input switched to continuous | Initial rate, localized swap | First mean, localized swap | First mean, full-region swap |
| --- | ---: | ---: | ---: |
| Parent | .07240916386 | .06601221489 | .06603842258 |
| Daughter | .06312562923 | .06693240436 | .06619076971 |
| Both | .05689879272 | .06600201819 | .06528878623 |

The native initial and first mean rates are .07863619025 and .06694243427. The local daughter change reproduces nearly all **instantaneous** daughter suppression but gives only 1.33% of the full daughter-swap reduction in the first mean. The local parent swap retains the full finite parent effect to within about 3% of the effect itself. Switching both only locally retains about 56.9% of the full finite two-species contrast, even though the starting rate is nearly the same. This establishes sensitivity to the broader starting hard distribution during this joined evolution, not a measured transport flux from any particular outer cell. All accepted histories pass the unchanged local and terminal energy and terminal entropy gates.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
for p in local-parent local-daughter local-both; do
  python -W error early_pole_prepared_hard_swap_trajectory_reference.py "$a" --case "$p" --steps 16
done
```

Raw: `results/early-pole-prepared-local-hard-swaps-2026-11-27.json.gz`.
