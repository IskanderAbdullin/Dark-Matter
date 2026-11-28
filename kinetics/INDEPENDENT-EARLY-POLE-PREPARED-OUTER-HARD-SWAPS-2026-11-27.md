# The finite daughter effect is carried by changes outside the initial cold aperture

Use the same native first-window model and reinitialize the active response at each substituted state. For each charged species, leave the [two cold-active cells](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-KINEMATIC-APERTURE-2026-11-27.md) native and substitute the continuous preparation in the other 18 retained hard cells. This complements the [localized swaps](INDEPENDENT-EARLY-POLE-PREPARED-LOCAL-HARD-SWAPS-2026-11-27.md); an outside moment change can still shift the common dispersion. Compare each with the original native state and with the whole-species swap at the same 16-step width.

| Continuous input outside the two-cell aperture | Initial cold rate | First mean rate | Fraction of full-species mean reduction |
| --- | ---: | ---: | ---: |
| None | .07863619025 | .06694243427 | — |
| Parent only | .07863951819 | .06696473823 | −2.47% |
| Daughter only | .07863573739 | .06620059799 | 98.69% |
| Both | .07863906534 | .06622334784 | 43.48% |

The outside-only daughter swap initially changes the rate by only −4.53e−7 but changes the first mean by −.00074184 with the hard channels joined. Turn off both elastic channels while keeping reciprocal decay, photon pairs, Compton and flavor conversion. Relative to that independently evolved nonelastic native reference, the same outside daughter input changes the first mean by only −4.71e−7, about 0.064% of the joined effect. The substantially larger finite-time response thus requires the elastic channels in this controlled model, although this does not isolate a unique microscopic transport path or exclude effects on the mean field. All trajectories pass the original gates; no statement is made about hard cells outside the retained box or physical abundance.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
p=early_pole_prepared_hard_swap_trajectory_reference.py
for k in outside-parent outside-daughter outside-both; do
  python -W error "$p" "$a" --case "$k" --steps 16
done
python -W error "$p" "$a" --case native --steps 16 --no-elastic
python -W error "$p" "$a" --case outside-daughter --steps 16 --no-elastic
```

Raw: `results/early-pole-prepared-outer-hard-swaps-2026-11-27.json.gz`.
