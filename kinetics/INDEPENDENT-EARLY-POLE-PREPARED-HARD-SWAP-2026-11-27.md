# The daughter reconstruction accounts for most of the instantaneous prepared-rate gap

The two prepared fine hard states have identical cell edges, identical photon state and identical initial cold scalar bins. Hold that support and the native archived exterior fixed; replace only the parent or only the daughter retained hard moments with the count-conserving continuous reconstruction. Replacing both reproduces every mode of the independently recorded continuous starting cold source. We evaluate the reciprocal decay source at those prepared states and do not evolve them in this check.

| Retained hard moments | Instantaneous cold rate |
| --- | ---: |
| Native parent and daughter | .07863619025 |
| Continuous parent only | .07241425886 |
| Continuous daughter only | .06312505721 |
| Both continuous | .05690395094 |
| Continuous parent retuned to native parent moments | .07283452112 |
| Retuned parent plus continuous daughter | .05732364429 |

The individual swaps account for 28.63% and 71.37% of the starting signed contrast, with a small static interaction of +8.25e−7 in the rate. The raw parent swap does not conserve the number of parents **within the retained hard window**: its count changes from .292439 to .294240, although the term “count-conserving” describes the original reconstruction. Retune the substituted parent using two shared affine changes to its FD logits to match the native retained parent number and *bare* energy. The residuals are at most 1.2e−16. This manufactured moment match leaves 98.1% of the original starting contrast when the continuous daughter is also substituted. Matching bare energy does not match all interaction moments or the full spectrum.

The preceding [first-window time study](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-TIME-LOCAL-MEMORY-2026-11-27.md) distinguishes the large instantaneous starting contrast from a much smaller mean-rate difference after relaxation. These swaps only attribute the starting contrast; they neither determine which spectrum is physical nor give a channel attribution integrated over that window.

From `kinetics/`:
```sh
OPENBLAS_NUM_THREADS=1 python -W error early_pole_prepared_hard_swap_reference.py \
  results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
```

Raw: `results/early-pole-prepared-hard-swap-2026-11-27.json`.
