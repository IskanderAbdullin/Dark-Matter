# Time convergence of a coupled, entropy-limited two-mass chart

Use the same [entropy-coordinate event chart](INDEPENDENT-TWO-MASS-ENTROPY-EXTENT-2026-11-27.md)
with three admissible event directions sharing parent and daughter
populations. Compare two split time updates against an unsplit
adaptive solution of the sum of the same three manufactured rates.

The first method composes backward event steps. The second uses
implicit midpoint event steps in symmetric ordering. For each
midpoint substep, require a solution between the current state
and the unique entropy stationary point along that event
extent. This keeps every substep inside its positive/Pauli
bounds and makes entropy nondecreasing. If the midpoint
root is not bracketed, fall back to the entropy-bounded
backward step instead of extrapolating. The fallback need not
be second order.

| Fixture | Observed p, backward split | Observed p, symmetric midpoint | Maximum count error at 64 steps, backward / symmetric |
|---|---:|---:|---:|
| Ordinary | .9880 | 1.9998 | .0015624 / .00001235 |
| Near rest | .9882 | 1.9998 | .0013699 / .00001592 |

The rates are estimated from the last two mesh levels, against an
independently tolerance-refined unsplit DOP853 solution. Neither
four-level convergence sweep triggered the fallback. A separate
stiff control starting near the Pauli and positivity limits does
trigger it and increases entropy. The number and split energy
ledgers close on each solved step. A central finite difference of
the full entropy agrees with the event affinity to 1.85e−12 in
the two directional checks.

This tests time integration on a fixed manufactured chart with
chosen capacities and rates. It is not a convergence proof for
arbitrary time steps, the actual physical evolution, or a prediction
of occupations or transfer functions.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_entropy_coupled_time_reference.py
```
[Full receipt](results/two-mass-entropy-coupled-time-2026-11-27.json).
