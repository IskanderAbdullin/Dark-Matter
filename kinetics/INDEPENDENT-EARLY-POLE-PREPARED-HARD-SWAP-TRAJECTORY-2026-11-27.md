# The finite cold-source contrast is split more evenly between the two hard inputs

Evolve the [parent and daughter input swaps](INDEPENDENT-EARLY-POLE-PREPARED-HARD-SWAP-2026-11-27.md) over the first coarse interval (3.3333333333333335e−7) with 16 macro steps. Within this finite hard model, build the active reference response from the substituted prepared state so it does not include a spurious pre-swap response. The charged cells, photon data, initial cold scalar bins, interaction strengths and solver budget are otherwise fixed. The both-swapped count history has the same count at every step as the separately initialized continuous branch.

| Prepared changes on a native hard support | Instantaneous cold rate at t=0 | First interval mean |
| --- | ---: | ---: |
| Neither | .07863619025 | .06694243427 |
| Parent only | .07241425886 | .06603842258 |
| Daughter only | .06312505721 | .06619076971 |
| Both | .05690395094 | .06528878623 |
| Parent moment-matched to native, daughter left native | .07283452112 | .06643447117 |

The daughter-only swap accounts for 71.4% of the instantaneous starting difference, but 45.5% of the finite mean-rate difference. The raw parent-only contribution becomes 54.7% over the finite window. Retuning its retained number and bare energy to their native values reduces the parent-only contribution to 30.7% of the original finite difference. This retuning is a manufactured two-moment change, not a unique decomposition into normalizations. The signed finite parent–daughter interaction in the mean is +2.03e−6, about .12% of the net contrast: the two isolated substitutions nearly add *as initial-condition changes*. This should not be conflated with the [strong non-additivity of the elastic channel toggles](INDEPENDENT-EARLY-POLE-FIRST16-ELASTIC-FACTORIAL-2026-11-27.md).

All accepted paths satisfy the same terminal energy, entropy and elastic budget checks. The first-interval mean still has a [documented time-width dependence](INDEPENDENT-EARLY-POLE-FIRST-WINDOW-TIME-LOCAL-MEMORY-2026-11-27.md). The result does not determine the physical production spectrum or the asymptotic contributions of the two charged species.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
for c in native parent daughter both matched-parent; do
  python -W error early_pole_prepared_hard_swap_trajectory_reference.py "$a" --case "$c" --steps 16
done
```

Raw: `results/early-pole-prepared-hard-swap-trajectory-2026-11-27.json.gz`.
