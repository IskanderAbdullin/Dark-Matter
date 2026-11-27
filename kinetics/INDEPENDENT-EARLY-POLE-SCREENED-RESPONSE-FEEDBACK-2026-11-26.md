# Retarded active-occupation feedback in the longer screened reference

The [implicit finite-packet solver](INDEPENDENT-EARLY-POLE-SCREENED-DISCRETE-GRADIENT-2026-11-26.md)
now resolves in this finite model a reference interval `2e-6`,
25 times the first fully matched `8e-8` screened experiment.
The clock, finite radial domain, scalar support and separate
archived exterior are unchanged. Photon pairs are off in the
following dynamic comparison. No expansion or new gauge channel
is inferred from this extension.

A tighter nonlinear residual ensures that the energy ledger
continues to close over the extended stiff evolution.
At 24 backward split steps the relative effect of turning on
the screened elastic channel is:

| Initial occupation representation | Fine archive | Third archive |
| --- | ---: | ---: |
| Native steps | −9.89% | +10.01% |
| Native count-conserving continuous | +20.98% | −7.95% |

Halving the step size from 12 to 24 changes the actual gain
by no more than 0.34% on these tests. It does not remove the
representation-dependent sign. All sampled gains continue to
increase, with finite variational energy and entropy closure.
The long calculation therefore does not identify a freeze-out.

## Keep the exterior frozen; update the active polarization

Each active Fermi cell gives an isotropic derivative contribution
from the continuous logit profile and a jump at either boundary.
For the **difference from the same initial projected state**,
both are propagated through the longitudinal and transverse
massive Vlasov kernels. The exterior is left untouched: the
response is the archived polarization plus just this active
increment. The boundary contribution can be spread over a
small fraction of the adjacent cell widths while preserving
its static strength exactly. Independently, integration by
parts on the volume quadrature gives the same static Debye shift
using the factor `(E + p²/E)/p²`. No artificial shift occurs at
the unchanged initial state.

At each accepted coupled macro step, use the current state
just before the elastic step. Below the same matching transfer
`Q=.08`, update the real and imaginary longitudinal and transverse
denominators and their interference, by the ratio of updated to
original integrated current traces. Above the cut retain the
vacuum weights. All events remain paired with their reverse.
The coefficients are held within the implicit step: this is
a lagged collisionless comparison, not an iteration to a
self-consistent quantum plasma solution. A second comparator
changes only the Debye amplitude while keeping its original
frequency shape.

The Debye mass **squared** has changed by about +0.16% at the
last refresh. The longitudinal/transverse increment as a function
of phase velocity can be much larger in a few tabulated components,
but the modified individual soft packet weights differ by less
than 0.35% at that stage. The cold scalar gain differs from
the frozen-polarization calculation by only −0.000174% (fine)
and −0.000054% (third), at matched 12-step exposure. A Debye-only
rescaling is almost indistinguishable for that observable.
This is much smaller than the measured step-size and initial
representation differences, so the numbers are a sensitivity
check rather than a precision physical effect.

A separate same-state reweighting check after an implicit
`4e-7` elastic step varied edge spreading from zero to .25 and
doubled the phase grid. The pointwise increment is sensitive
to this choice, but the conditional integrated daughter heat
changes by less than `5e-7` between these updated kernels.
The full update, relative to holding the field response fixed,
changes that heat by about .0094 (against a frozen-kernel
heat of 257 on that test).

No conclusion about the final scalar distribution, a complete
cosmological rate, or a transfer function follows. The next
constraint is the narrow charged input distribution and the
rest of the gauge collision channels, not merely refining the
feedback described here.

From `kinetics/` with one BLAS thread:

```sh
python early_pole_screened_active_response_update_reference.py fine --duration 4e-7 --smear .25 --resolution 2
python early_pole_cold_screened_split_reference.py fine --duration 2e-6 --count 12 --screening-update full
python early_pole_cold_screened_split_reference.py fine --duration 2e-6 --count 24 --projection cell-conservative
```

Raw: `results/early-pole-screened-response-feedback-2026-11-26.json.gz`.
