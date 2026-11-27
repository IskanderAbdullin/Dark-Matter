# The daughter occupation contribution reverses during the first cold-rate relaxation

Along the [accepted matched 64-step first window](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-TIME-LOCAL-MEMORY-2026-11-27.md), replay each saved snapshot with the same cold packet integral. At each time, start from the native fitted moments, then substitute the continuous parent and daughter occupations at held native dispersion and scalar seed. Finally switch the hard energies and then the scalar seed to the continuous values. The first and last states of this ordered decomposition reproduce the saved cold-mode rates. These substitutions are diagnostics on recorded states, not counterfactual trajectories.

| Exposure | Parent occupation contribution | Daughter occupation contribution | Actual continuous minus native rate |
| ---: | ---: | ---: | ---: |
| 0 | −.00622702 | −.01551119 | −.02173224 |
| 4.1667e−8 | −.00140792 | −.00032146 | −.00172351 |
| 8.3333e−8 | −.00084857 | +.00013104 | −.00071285 |
| 3.3333e−7 | −.00019297 | +.00006884 | −.00012175 |

The initially dominant daughter contribution decreases quickly and reverses sign between the two listed samples. At the end it offsets some of the remaining negative parent contribution. Switching the energies at the last sample contributes +2.37e−6; switching the evolved cold seed then contributes only +7.36e−9. The two original contributing occupation cells per species, spanning approximately [1.4139, 1.4667], still account for all of the evaluated occupation-induced contrast at these nine snapshots. This locality concerns the instantaneous cold packet support, not where collisions transport hard occupations during the run. The finite orientation must not be extrapolated to later times or to a physical production spectrum.

From `kinetics/`:
```sh
OPENBLAS_NUM_THREADS=1 python -W error early_pole_first_window_occupation_history_reference.py \
  results/early-pole-first-window-64-memory-2026-11-27.json.gz
```

Raw: `results/early-pole-first-window-occupation-history-2026-11-27.json`.
